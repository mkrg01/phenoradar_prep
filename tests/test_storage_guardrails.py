"""Real child processes, flock leases, and verified maintenance protect storage."""
import gzip
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import fcntl

import pytest
import yaml

from common import write_json
import scratch_reservation as admission
from storage_runtime import settings, task_settings, storage_policy, environment
from storage_meter import monitored, append_sample, compact_series
from derived_cache_cleanup import retire
from test_storage_meter import settings as monitor_settings
from test_datasets import dataset_project, new_dataset, native_events


def test_admission_serializes_processes_and_reaps_dead_leases(tmp_path, monkeypatch):
    directory = tmp_path / 'leases'; directory.mkdir()
    first, second = tmp_path / 'first', tmp_path / 'second'
    first.mkdir(); second.mkdir()
    marker = first / 'restart-data'; marker.write_text('retain after crash')
    monkeypatch.setattr(admission, 'coordination_directory', lambda device: directory)
    monkeypatch.setattr(admission, 'free_space', lambda path: {'bytes': 1000000000, 'inodes': 1000})
    config = settings({'scratch_reserve_bytes': 600000000, 'scratch_reserve_inodes': 600})
    child = '''import json,sys,time
from pathlib import Path
import scratch_reservation as a
a.coordination_directory=lambda device:Path(sys.argv[1])
a.free_space=lambda path:{'bytes':1000000000,'inodes':1000}
with a.reserve(sys.argv[2],json.loads(sys.argv[3])):
 print('reserved',flush=True)
 time.sleep(30)
'''
    env = dict(os.environ, PYTHONPATH=str(Path(admission.__file__).parent))
    process = subprocess.Popen([sys.executable, '-c', child, str(directory), str(first), json.dumps(config)],
                               env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    try:
        assert process.stdout.readline().strip() == 'reserved'
        with pytest.raises(admission.ScratchAdmissionError, match='reservation unavailable'):
            with admission.reserve(second, config):
                pytest.fail('overcommitted')
        process.kill(); process.wait(timeout=5)
        with admission.reserve(second, config):
            assert len(list(directory.glob('lease-*.json'))) == 1
        assert not list(directory.glob('lease-*.json'))
        assert marker.read_text() == 'retain after crash'
    finally:
        if process.poll() is None:
            process.kill(); process.wait(timeout=5)


@pytest.mark.parametrize('config', [
    {'scratch_max_inodes': 0}, {'scratch_reserve_bytes': True},
    {'start_fraction': .98, 'stop_fraction': .95},
    {'scratch_reserve_bytes': 10, 'scratch_max_bytes': 20},
    {'task_limits': {'sample': {'max_bytes': 1}}},
])
def test_invalid_limits_are_rejected(config):
    with pytest.raises(ValueError):
        settings(config)


def test_full_scratch_filesystem_blocks_before_creating_work(tmp_path, monkeypatch):
    import storage_runtime
    root, node = tmp_path / 'project', tmp_path / 'node'
    root.mkdir(); node.mkdir()
    def full(*args, **kwargs):
        raise OSError(28, 'No space left on device')
    monkeypatch.setattr(storage_runtime.tempfile, 'mkdtemp', full)
    with pytest.raises(admission.ScratchAdmissionError, match='No space left'):
        with storage_runtime.task_scratch(root, settings({'scratch_root':str(node)})):
            pytest.fail('child could start on full scratch')


def test_both_scheduler_command_styles_select_the_same_task_policy():
    storage = {'task_limits': {'odb_map': {'scratch_reserve_bytes': 123456, 'scratch_max_inodes': 500}}}
    age = task_settings(storage, 'odb_map', ['bash', 'jobscript.sh'])
    slurm = task_settings(storage, 'workflow', ['bash', '-c', 'python -m snakemake --target-jobs odb_map:chunk=chunk_000 --cores 2'])
    assert age == slurm
    assert age['scratch_max_bytes'] == age['scratch_reserve_bytes'] == 123456


def test_retry_policy_is_storage_only_and_records_content_identity(tmp_path, monkeypatch):
    path = tmp_path / 'policy.yaml'
    path.write_text('storage:\n  task_limits:\n    sample:\n      scratch_max_bytes: 2000\n      scratch_reserve_bytes: 2500\n')
    monkeypatch.setenv('PHENORADAR_STORAGE_POLICY', str(path))
    original = {'task_limits': {'sample': {'scratch_min_free_inodes': 100}}}
    merged, receipt = storage_policy(original, tmp_path)
    assert receipt['sha256'] and merged['task_limits']['sample']['scratch_max_bytes'] == 2000
    assert merged['task_limits']['sample']['scratch_min_free_inodes'] == 100
    assert 'scratch_max_bytes' not in original['task_limits']['sample']
    path.write_text('storage: {}\ntranslation: {table: 2}\n')
    with pytest.raises(ValueError, match='only a storage'):
        storage_policy(original, tmp_path)


@pytest.mark.parametrize('variables,threads', [({'NSLOTS': '3'}, '3'),
    ({'PHENORADAR_TASK_CPUS': '2', 'NSLOTS': '8'}, '2'),
    ({'SLURM_CPUS_PER_TASK': '4', 'PHENORADAR_TASK_CPUS': '16'}, '4')])
def test_blas_threads_follow_allocations_inside_cleanenv(tmp_path, monkeypatch, variables, threads):
    for key in ('NSLOTS', 'PHENORADAR_TASK_CPUS', 'SLURM_CPUS_PER_TASK'):
        monkeypatch.delenv(key, raising=False)
    for key, value in variables.items():
        monkeypatch.setenv(key, value)
    env = environment(tmp_path, settings({}))
    assert env['OPENBLAS_NUM_THREADS'] == threads
    assert env['APPTAINERENV_OMP_NUM_THREADS'] == threads
    assert env['SINGULARITYENV_MKL_NUM_THREADS'] == threads


def test_cli_capacity_denial_never_starts_the_child(tmp_path):
    root, node = tmp_path / 'project', tmp_path / 'node'
    root.mkdir(); node.mkdir()
    config = root / 'pipeline.yaml'
    config.write_text(yaml.safe_dump({'storage': {'scratch_root': str(node), 'scratch_reserve_bytes': 10**30}}))
    marker = root / 'child-started'
    script = Path(__file__).parents[1] / 'workflow/scripts/storage_runtime.py'
    result = subprocess.run([sys.executable, str(script), 'run', '--pipeline', str(config), '--root', str(root),
                             '--usage', str(root / 'jobs/storage'), '--label', 'sample', '--', sys.executable,
                             '-c', 'from pathlib import Path; Path(' + repr(str(marker)) + ').touch()'],
                            capture_output=True, text=True, timeout=10)
    assert result.returncode and not marker.exists()
    receipt = json.loads(next((root / 'jobs/storage').glob('*.json')).read_text())
    assert receipt['state'] == 'blocked' and receipt['child_started'] is False
    assert not list(node.iterdir())


def test_scratch_inode_limit_stops_a_real_writer(tmp_path):
    root, scratch = tmp_path / 'project', tmp_path / 'scratch'
    root.mkdir(); scratch.mkdir()
    receipt = root / 'usage.json'
    command = [sys.executable, '-c', 'import pathlib,sys,time; p=pathlib.Path(sys.argv[1]); '
               '[ (p/str(i)).touch() for i in range(100) ]; time.sleep(30)', str(scratch)]
    with pytest.raises(ValueError, match='scratch inodes'):
        monitored(command, cwd=root, env=dict(os.environ), shared_root=root, scratch=scratch,
                  sample_roots=[], receipt=receipt, settings=dict(monitor_settings(), scratch_max_inodes=20))
    saved = json.loads(receipt.read_text())
    assert saved['state'] == 'failed' and saved['child_exit_code'] != 0
    assert saved['scratch_peak']['inodes'] >= 20


def test_shared_watchdog_stops_with_headroom_before_the_absolute_budget(tmp_path):
    root, scratch = tmp_path / 'project', tmp_path / 'scratch'
    root.mkdir(); scratch.mkdir(); receipt = root / 'usage.json'
    command = [sys.executable, '-c', 'import pathlib,sys,time; pathlib.Path(sys.argv[1]).write_bytes(b"x"*1900000); time.sleep(30)', str(root / 'data')]
    with pytest.raises(ValueError, match='budget reached'):
        monitored(command, cwd=root, env=dict(os.environ), shared_root=root, scratch=scratch,
                  sample_roots=[], receipt=receipt, settings=dict(monitor_settings(), max_bytes=2000000, stop_fraction=.95))
    peak = json.loads(receipt.read_text())['shared_peak']['allocated_bytes']
    assert 1900000 <= peak < 2000000


def test_background_writer_does_not_outlive_its_owner(tmp_path):
    root, scratch = tmp_path / 'project', tmp_path / 'scratch'
    root.mkdir(); scratch.mkdir(); marker = scratch / 'writer'
    writer = 'import pathlib,time; p=pathlib.Path(' + repr(str(marker)) + ');\nwhile True:\n p.write_text(str(time.time())); time.sleep(.01)'
    command = [sys.executable, '-c', 'import subprocess,sys,time; subprocess.Popen([sys.executable,"-c",sys.argv[1]]); time.sleep(.2)', writer]
    monitored(command, cwd=root, env=dict(os.environ), shared_root=root, scratch=scratch,
              sample_roots=[], receipt=root / 'usage.json', settings=monitor_settings())
    before = marker.read_bytes(); time.sleep(.1)
    assert marker.read_bytes() == before


def test_time_series_has_bounded_storage_and_compresses_at_completion(tmp_path):
    receipt = tmp_path / 'usage.json'
    for i in range(1000):
        append_sample(receipt, {'at': str(i), 'data': 'number=' + str(i) * 100}, 1024)
    assert receipt.with_suffix('.jsonl').stat().st_size <= 1024
    assert len(list(tmp_path.iterdir())) <= 2
    compact_series(receipt)
    assert not receipt.with_suffix('.jsonl').exists()
    assert gzip.decompress(receipt.with_suffix('.jsonl.gz').read_bytes()).strip()
    assert sum(p.stat().st_size for p in tmp_path.iterdir()) < 2048


def test_private_derivations_retire_but_published_snapshot_links_survive(tmp_path):
    cache = tmp_path / 'work/cache/odb'
    digest = 'a' * 64
    target = cache / '.tables' / digest
    target.mkdir(parents=True); table = target / 'table.tsv.gz'; table.write_bytes(b'published table')
    snapshot = cache / 'v12_3193/run'; snapshot.mkdir(parents=True)
    os.link(table, snapshot / 'table.tsv.gz')
    write_json(snapshot / 'snapshot.json', {'table': {'path': 'table.tsv.gz'}})
    lock = target.with_suffix('.lock'); lock.touch()
    with lock.open('r') as handle:
        fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        with pytest.raises(BlockingIOError):
            retire(tmp_path, [], apply=True)
    assert retire(tmp_path, [], apply=False)['targets']
    retire(tmp_path, [], apply=True)
    assert not target.exists() and lock.exists()
    assert (snapshot / 'table.tsv.gz').read_bytes() == b'published table'


def test_successful_worker_compacts_before_whole_build_completes(dataset_project):
    from dataset import load, stage_workspace, worker
    build = new_dataset(dataset_project)
    stage_workspace(build, load(build), [1])
    work = build / 'work/genegalleon/New_plant_SRR1'
    cache = work / '.gg_cache/content_digests.sqlite3'
    cache.parent.mkdir(); cache.write_bytes(b'cached digest row\n' * 10000)
    worker(build, 1)
    assert not (build / 'completed.json').exists()
    assert not cache.exists()
    assert (work / 'diagnostics.tar.gz').is_file()
    events = native_events(build)
    worker(build, 1)
    assert native_events(build) == events
    saved = json.loads((build / 'jobs/status/New_plant_SRR1.sample.json').read_text())
    assert saved['native_exit_code'] == 0
