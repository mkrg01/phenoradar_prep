"""Native evidence proves checkpoints; scheduler failures remain separate evidence."""
import fcntl
import json
import os
import subprocess
from pathlib import Path

import pytest

from dataset import cleanup_genegalleon, load, status, submit, worker
from genegalleon_run import observation, scheduler_states
from storage_management import run_storage
from test_datasets import dataset_project, new_dataset, native_events


def test_busco_failure_restarts_busco_without_assembling_again(dataset_project, monkeypatch):
    build = new_dataset(dataset_project)
    submit(build, until='quant', dry_run=True)
    monkeypatch.setenv('FAKE_GG_FAIL_BUSCO', '1')
    with pytest.raises(subprocess.CalledProcessError): worker(build, 1)
    row = status(build)[0]
    assert [row[s] for s in ('assembly', 'busco', 'quant')] == ['reuse', 'pending', 'pending']
    assert row['stopped_at'] == 'busco'
    monkeypatch.delenv('FAKE_GG_FAIL_BUSCO')
    worker(build, 1)
    assert native_events(build)[-1]['env']['GG_TRANSCRIPTOME_RUN_ASSEMBLY'] == '0'
    assert status(build)[0]['stopped_at'] is None


@pytest.mark.parametrize('response', ['{"schema":"unsupported"}', '[]'])
def test_unavailable_native_api_never_proves_partial_completion(dataset_project, monkeypatch, response):
    build = new_dataset(dataset_project)
    submit(build, until='quant', dry_run=True)
    run = subprocess.run
    def unsupported(command, **kwargs):
        if any(str(arg).endswith('workflow_api.py') for arg in command):
            return subprocess.CompletedProcess(command, 0, response, '')
        return run(command, **kwargs)
    monkeypatch.setattr(subprocess, 'run', unsupported)
    monkeypatch.setenv('FAKE_GG_FAIL_BUSCO', '1')
    with pytest.raises(subprocess.CalledProcessError): worker(build, 1)
    assert status(build)[0]['assembly'] == 'pending'
    monkeypatch.delenv('FAKE_GG_FAIL_BUSCO')
    worker(build, 1)
    steps = [json.loads(line) for line in (build / 'work/genegalleon/New_plant_SRR1/steps.jsonl').read_text().splitlines()]
    # The native engine reused its own checkpoint, even without a wrapper receipt.
    assert steps.count('transcriptome_longest_cds') == 1
    assert status(build)[0]['quant'] == 'reuse'


def test_cleanup_retains_scratch_held_by_orphan_native_process(dataset_project):
    build = new_dataset(dataset_project)
    submit(build, until='quant', dry_run=True)
    manifest = load(build)
    scratch = build / 'work/genegalleon/New_plant_SRR1/output/transcriptome_assembly/tmp/1_native'
    scratch.mkdir(parents=True)
    (scratch / 'large').write_text('still in use')
    with (scratch / '.gg_active.lock').open('w') as handle:
        fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        with pytest.warns(RuntimeWarning, match='cleanup incomplete'):
            result = cleanup_genegalleon(build, manifest, manifest['items'][0], 'sample', {})
        assert result['state'] == 'pending' and (scratch / 'large').exists()
    result = cleanup_genegalleon(build, manifest, manifest['items'][0], 'sample', {})
    assert result['state'] == 'complete' and not scratch.exists()


def test_accounting_queries_exact_array_tasks_and_ignores_steps(monkeypatch):
    calls = []
    def accounting(command, **kwargs):
        calls.append(command)
        return subprocess.CompletedProcess(command, 0, '1001_2|OUT_OF_MEMORY|0:9\n1001_2.batch|FAILED|0:9\n1001_2||\n', '')
    monkeypatch.setattr(subprocess, 'run', accounting)
    states = scheduler_states([{'slurm_id': '1001_2', 'job_id': '9999'}, {'job_id': None}])
    assert states == {'1001_2': {'state': 'OUT_OF_MEMORY', 'exit_code': '0:9'}}
    assert calls[0][calls[0].index('--jobs') + 1] == '1001_2'


def test_status_uses_scheduler_evidence_after_wrapper_oom(dataset_project, monkeypatch):
    build = new_dataset(dataset_project)
    submit(build, until='quant', dry_run=True)
    monkeypatch.setenv('SLURM_JOB_ID', '1001')
    monkeypatch.setenv('SLURM_ARRAY_JOB_ID', '1000')
    monkeypatch.setenv('SLURM_ARRAY_TASK_ID', '2')
    monkeypatch.setenv('FAKE_GG_FAIL_QUANT', '1')
    with pytest.raises(subprocess.CalledProcessError): worker(build, 1)
    receipt = build / 'jobs/status/New_plant_SRR1.sample.json'
    job = json.loads(receipt.read_text())
    job['state'] = 'running'  # SIGKILL prevented the final wrapper record.
    receipt.write_text(json.dumps(job))
    run = subprocess.run
    def accounting(command, **kwargs):
        if command[0] == 'sacct':
            return subprocess.CompletedProcess(command, 0, '1000_2|OUT_OF_MEMORY|0:9\n', '')
        return run(command, **kwargs)
    monkeypatch.setattr(subprocess, 'run', accounting)
    row = status(build)[0]
    assert row['jobs']['sample']['slurm_id'] == '1000_2'
    assert row['jobs']['sample']['state'] == 'failed'
    assert row['jobs']['sample']['scheduler']['state'] == 'OUT_OF_MEMORY'
    assert row['stopped_at'] == 'quant'


@pytest.mark.parametrize('state', ['failed', 'running'])
def test_inactive_failed_sample_scratch_can_be_cleaned_without_computation(dataset_project, monkeypatch, state):
    build = new_dataset(dataset_project)
    submit(build, until='quant', dry_run=True)
    monkeypatch.setenv('FAKE_GG_FAIL_QUANT', '1')
    with pytest.raises(subprocess.CalledProcessError): worker(build, 1)
    receipt = build / 'jobs/status/New_plant_SRR1.sample.json'
    job = json.loads(receipt.read_text())
    job['state'] = state
    receipt.write_text(json.dumps(job))
    work = build / 'work/genegalleon/New_plant_SRR1/output/transcriptome_assembly'
    scratch = work / 'tmp/1_native'
    scratch.mkdir(parents=True)
    (scratch / 'large').write_bytes(b'x' * 8192)
    reads = work / 'amalgkit_getfastq/New_plant/SRR1/reads.fastq.gz'
    reads.parent.mkdir(parents=True)
    reads.write_bytes(b'completed download for retry')
    before = native_events(build)
    preview = run_storage(build, 'build', inspect=True)
    assert preview['counts'] == {'available': 1}
    assert preview['jobs'][0]['job_state'] == state
    assert scratch.exists() and reads.exists()
    manifest = load(build)
    with (Path(manifest['config']['store']) / 'New_plant_SRR1/.worker.lock').open('r') as handle:
        fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        assert run_storage(build, 'build', inspect=True, apply=True)['counts'] == {'blocked': 1}
        assert scratch.exists()
    with (scratch / '.gg_active.lock').open('w') as handle:
        fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        assert run_storage(build, 'build', inspect=True, apply=True)['counts'] == {'blocked': 1}
        assert scratch.exists()
    assert run_storage(build, 'build', inspect=True, apply=True)['counts'] == {'complete': 1}
    assert not scratch.exists() and reads.exists()
    assert native_events(build) == before


@pytest.mark.parametrize('stage,function', [('assembly', 'register_reference'), ('busco', 'register_busco')])
def test_publication_failure_reports_the_unpublished_stage(dataset_project, monkeypatch, stage, function):
    import dataset
    build = new_dataset(dataset_project)
    submit(build, until='quant', dry_run=True)
    register = getattr(dataset, function)
    def fail(*args, **kwargs):
        raise ValueError('publication failed')
    monkeypatch.setattr(dataset, function, fail)
    with pytest.raises(ValueError, match='publication failed'): worker(build, 1)
    row = status(build)[0]
    assert row['jobs']['sample']['phase'] == 'publication'
    assert row['stopped_at'] == stage
    monkeypatch.setattr(dataset, function, register)
    worker(build, 1)
    assert status(build)[0]['quant'] == 'reuse'
    assert len(native_events(build)) == 1


def test_malformed_observation_is_unavailable(tmp_path):
    attempt = tmp_path / 'output/observations/bad'
    attempt.mkdir(parents=True)
    (attempt / 'run.json').write_text('[]')
    assert observation(tmp_path)['state'] == 'unavailable'


def test_worker_passes_the_native_lock_to_its_launcher(dataset_project, monkeypatch):
    build = new_dataset(dataset_project)
    submit(build, until='quant', dry_run=True)
    run = subprocess.run
    inherited = []
    def launch(command, **kwargs):
        if command[0] == 'bash' and str(command[1]).endswith('gg_transcriptome_generation_entrypoint.sh'):
            descriptor, = kwargs['pass_fds']
            lock = build / 'work/genegalleon/New_plant_SRR1/.native.lock'
            assert os.fstat(descriptor).st_ino == lock.stat().st_ino
            inherited.append(descriptor)
        return run(command, **kwargs)
    monkeypatch.setattr(subprocess, 'run', launch)
    worker(build, 1)
    assert len(inherited) == 1
