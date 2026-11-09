"""Pilot selection and storage policies are shared by Slurm and AGE builds."""
import csv
import io
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest
import yaml

import dataset
import storage_runtime
from common import read_tsv, write_json, write_tsv
from test_datasets import dataset_project

ROOT = Path(__file__).resolve().parents[1]


def test_pilot_freezes_only_subset_without_a_site_specific_configuration(dataset_project, monkeypatch):
    root = dataset_project
    rows = [{'scientific_name': f'Plant species{i}', 'run': f'TEST{i}', 'taxid': str(i + 100),
             'total_bases': str(i * 10000)} for i in range(40)]
    source = root / 'input/metadata.tsv'
    write_tsv(source, list(rows[0]), rows)
    monkeypatch.setattr(dataset, 'resolve_software', lambda *a, **k: pytest.fail('submit host downloaded software'))
    build = dataset.prepare(root, None, root / 'config/build.yaml', defer=True, pilot_size=20)
    manifest = dataset.load(build)
    assert len(manifest['items']) == len(read_tsv(build / 'metadata.tsv')) == 20
    assert len(read_tsv(source)) == 40 and len(read_tsv(build / 'source_metadata.tsv')) == 40
    assert manifest['pilot']['source_samples'] == 40
    assert manifest['dataset_name'].endswith('_pilot')
    assert manifest['config']['reuse_from'] == []
    assert set(item['row']['run'] for item in manifest['items']) >= {'TEST0', 'TEST39'}
    assert not list(root.glob('config/*hgc*'))


def test_storage_preserves_reference_mounts_and_replaces_parent_temporary_path(tmp_path, monkeypatch):
    project, node = tmp_path / 'project', tmp_path / 'node'
    project.mkdir(); node.mkdir()
    config = storage_runtime.settings({'scratch_root': str(node)})
    monkeypatch.setenv('PHENORADAR_TASK_SCRATCH', '/old/node/task')
    monkeypatch.setenv('APPTAINER_MOUNT', 'type=bind,src=/ref,dst=/db,ro\ntype=bind,src=/old/node/task,dst=/old/node/task')
    with storage_runtime.task_scratch(project, config) as local:
        env = storage_runtime.environment(local, config)
        mounts = list(csv.reader(io.StringIO(env['APPTAINER_MOUNT'])))
        assert {len(row) for row in mounts} == {4}  # Apptainer's CSV parser requires equal widths.
        assert ['type=bind', 'src=/ref', 'dst=/db', 'ro'] in mounts
        assert not any('src=/old/node/task' in mount for mount in mounts)
        assert env['APPTAINERENV_PHENORADAR_TASK_SCRATCH'] == str(local)
        assert env['SINGULARITYENV_PHENORADAR_COMPACT_ODB_LOGS'] == '1'
    assert not local.exists() and node.is_dir()


def test_age_sparse_retry_array_preserves_sample_identity_and_measures_native_job(dataset_project, monkeypatch):
    import scheduler
    from test_datasets import new_dataset, native_events
    build = new_dataset(dataset_project, ('New plant a', 'New plant b', 'New plant c'))
    dataset.stage_workspace(build, dataset.load(build), [2])
    dataset.worker(build, 2)
    monkeypatch.setenv('PHENORADAR_SCHEDULER', 'age')
    monkeypatch.setattr(scheduler, 'age_capabilities', lambda: ('def_slot', 's_vmem'))
    monkeypatch.setattr(scheduler, 'array_size', lambda settings: 1000)
    commands = dataset.submit(build, until='quant', dry_run=True)
    assert len(commands) == 1 and commands[0][0] == 'qsub'
    assert commands[0][commands[0].index('-t') + 1] == '1-2'
    assert json.loads((build / 'jobs/0001_sample.indices.json').read_text()) == [1, 3]
    env = dict(os.environ, JOB_ID='88', SGE_TASK_ID='2', NSLOTS='4',
               PHENORADAR_TASK_CPUS='4', PHENORADAR_TASK_MEM_MB='128000')
    result = subprocess.run(['bash', commands[0][-1]], env=env, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    assert len(native_events(build)) == 2
    final = dataset.load(build)['items'][2]['species']
    receipt = json.loads((build / 'jobs/status' / (final + '.sample.json')).read_text())
    assert receipt['scheduler_id'] == '88_2' and receipt['scheduler'] == 'age'
    metrics = json.loads(next((build / 'jobs/storage').glob('sample_*.json')).read_text())
    assert metrics['state'] == 'complete' and metrics['sample_peak']['inodes'] > 0
    (build / 'jobs/0001_sample.indices.json').write_text('[3,1]')
    rejected = subprocess.run(['bash', commands[0][-1]], env=env, capture_output=True, text=True)
    assert rejected.returncode != 0 and len(native_events(build)) == 2


def test_common_storage_cli_measures_a_real_child_and_cleans_only_owned_scratch(tmp_path):
    project = tmp_path / 'project'
    project.mkdir()
    node = tmp_path / 'node'
    node.mkdir()
    (node / 'unrelated').write_text('keep')
    pipeline = project / 'pipeline.yaml'
    pipeline.write_text(yaml.safe_dump({'storage': {'scratch_root': str(node), 'interval_seconds': 1}}))
    command = [sys.executable, str(ROOT / 'workflow/scripts/storage_runtime.py'), 'run',
               '--pipeline', str(pipeline), '--root', str(project), '--usage', str(project / 'jobs/storage'),
               '--label', 'sample', '--', sys.executable, '-c',
               'import os,pathlib; pathlib.Path(os.environ["PHENORADAR_TASK_SCRATCH"],"payload").write_bytes(b"x"*10000)']
    subprocess.run(command, check=True)
    result = json.loads(next((project / 'jobs/storage').glob('*.json')).read_text())
    assert result['state'] == 'complete' and result['scratch_peak']['logical_bytes'] >= 10000
    assert list(node.iterdir()) == [node / 'unrelated']


def test_unfinished_storage_report_does_not_extrapolate(tmp_path, monkeypatch):
    build = tmp_path / 'results/pilot'
    build.mkdir(parents=True)
    monkeypatch.setattr(dataset, 'load', lambda *a, **k: {'root': str(tmp_path), 'config': {'storage': {}}, 'items': [{}] * 20})
    result = storage_runtime.report(build, 8000)
    assert result['projection'] is None and not result['complete']


def test_storage_projection_includes_native_work_and_checks_completed_database(tmp_path, monkeypatch):
    import build_products
    build = tmp_path / 'results/pilot'
    build.mkdir(parents=True)
    write_json(build / 'completed.json', {})
    (build / 'native_hits').write_bytes(b'x' * 100000)
    monkeypatch.setattr(dataset, 'load', lambda *a, **k: {'root': str(tmp_path), 'config': {'storage': {'max_bytes': 2500000000000}}, 'items': [{}] * 20})
    checked = []
    monkeypatch.setattr(build_products, 'load_complete', lambda path: checked.append(path))
    result = storage_runtime.report(build, 8000)
    assert checked == [build]
    assert result['retained_build']['allocated_bytes'] >= 100000
    assert result['projection']['allocated_bytes']['estimate'] >= result['retained_build']['allocated_bytes'] * 400
    assert result['projection']['allocated_bytes']['within_budget']
