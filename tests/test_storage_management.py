"""Batch cleanup verifies products, respects locks and never starts computation."""
import fcntl
import json
import subprocess
import sys
from pathlib import Path

import pytest

from common import write_json
from dataset import load, submit, worker
from run_odb_chunk import run as run_odb
from storage_management import run_storage, storage_report
from test_storage_cleanup import odb_job
from test_datasets import dataset_project, new_dataset, native_events

ROOT = Path(__file__).resolve().parents[1]


def test_batch_cleanup_preserves_failed_changed_and_active_jobs(odb_job, tmp_path):
    args, events = odb_job
    args['work_dir'] = tmp_path / 'work/orthogroups/mapping'
    pipeline = {'run_name': 'test', 'work_root': 'work', 'output_root': 'output'}
    manifest = {'root': str(tmp_path), 'config': {}, 'pipeline': pipeline}
    work = {}
    for n in range(3):
        label = f'chunk_{n:03d}'
        args.update(label=label, output_dir=tmp_path / 'output/orthogroups/mapping/chunks' / label)
        run_odb(**args, defer_cleanup=True)
        base = args['work_dir'] / label
        state = json.loads((base / 'status.json').read_text())
        work[label] = Path(state['work'])
        if n == 0:
            write_json(base / 'cleanup.json', {'state': 'pending', 'errors': [{'path': work[label].name}]})
        if n == 1:
            (args['output_dir'] / f'{label}.og.hits').write_text('changed')
        if n == 2:
            write_json(base / 'status.json', dict(state, state='failed'))
    before = events.read_text()
    preview = storage_report(tmp_path, manifest, 'analysis', inspect=True)
    assert preview['counts'] == {'pending': 1, 'blocked': 1, 'failed': 1}
    assert all(path.is_dir() for path in work.values())
    lock_path = args['work_dir'] / 'chunk_000/.lock'
    with lock_path.open('r') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        active = storage_report(tmp_path, manifest, 'analysis', inspect=True, apply=True)
        assert active['counts']['blocked'] == 2
        assert work['chunk_000'].is_dir()
    applied = storage_report(tmp_path, manifest, 'analysis', inspect=True, apply=True)
    assert applied['counts'] == {'complete': 1, 'blocked': 1, 'failed': 1}
    assert not work['chunk_000'].exists()
    assert work['chunk_001'].is_dir() and work['chunk_002'].is_dir()
    assert events.read_text() == before


def test_batch_cleanup_respects_debug_retention(odb_job, tmp_path):
    args, _ = odb_job
    args.update(work_dir=tmp_path / 'work/orthogroups/mapping',
                output_dir=tmp_path / 'output/orthogroups/mapping/chunks/chunk_000')
    run_odb(**args, defer_cleanup=True)
    manifest = {'root': str(tmp_path), 'config': {'storage': {'keep_intermediates': True}},
                'pipeline': {'run_name': 'test', 'work_root': 'work', 'output_root': 'output'}}
    report = storage_report(tmp_path, manifest, 'analysis', inspect=True, apply=True)
    assert report['counts'] == {'retained': 1}
    assert not (args['work_dir'] / 'chunk_000/cleanup.json').exists()


def test_build_status_and_cleanup_retry_without_rerunning_genegalleon(dataset_project):
    build = new_dataset(dataset_project)
    submit(build, until='quant', dry_run=True)
    for stage in ('assembly', 'busco', 'quant'):
        worker(build, stage, 1)
    before = native_events(build)
    work = build / 'work/genegalleon/New_plant_SRR1/output/transcriptome_assembly/tmp'
    work.mkdir(parents=True, exist_ok=True)
    (work / 'large').write_bytes(b'x' * 8192)
    receipt = build / 'jobs/cleanup/New_plant_SRR1.quant.json'
    write_json(receipt, {'state': 'pending', 'errors': [{'path': str(work), 'error': 'temporary'}]})
    status = subprocess.run([sys.executable, str(ROOT / 'workflow/scripts/dataset.py'), 'status',
                             '--build', str(build), '--storage'], text=True, capture_output=True)
    assert status.returncode == 0, status.stderr
    assert json.loads(status.stdout)['storage']['counts']['pending'] == 1
    before_receipt = receipt.read_bytes()
    preview = run_storage(build, 'build', inspect=True)
    pending = next(job for job in preview['jobs'] if job['state'] == 'pending')
    assert pending['inspection']['reclaimable_bytes'] >= 8192
    assert work.is_dir() and receipt.read_bytes() == before_receipt
    manifest = load(build)
    lock_path = Path(manifest['config']['store']) / 'New_plant_SRR1/.worker.lock'
    with lock_path.open('r') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        assert run_storage(build, 'build', inspect=True, apply=True)['counts']['blocked'] == 1
        assert work.exists()
    applied = subprocess.run([sys.executable, str(ROOT / 'workflow/scripts/dataset.py'), 'cleanup',
                              '--build', str(build), '--apply'], text=True, capture_output=True)
    assert applied.returncode == 0, applied.stderr
    assert not work.exists()
    assert json.loads(receipt.read_text())['state'] == 'complete'
    assert native_events(build) == before


def test_bulk_cleanup_never_follows_linked_job_parent(tmp_path):
    external = tmp_path / 'external'
    external.mkdir()
    (tmp_path / 'work').symlink_to(external, target_is_directory=True)
    manifest = {'root': str(tmp_path), 'config': {},
                'pipeline': {'run_name': 'test', 'work_root': 'work', 'output_root': 'output'}}
    with pytest.raises(ValueError, match='unsafe cleanup path'):
        storage_report(tmp_path, manifest, 'analysis', inspect=True, apply=True)
