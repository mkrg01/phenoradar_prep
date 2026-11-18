"""Stream completed worker logs into ordinary tar.gz archives without data loss."""
from concurrent.futures import ThreadPoolExecutor
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile

import pytest

import age_logs
from common import write_json

ROOT = Path(__file__).resolve().parents[1]
COMPLETE = {'state': 'COMPLETED', 'source': 'qacct', 'exit_code': '0:0'}


def worker(tmp_path, identifier='101'):
    archive = tmp_path / 'jobs/logs/age-workers-99.tar.gz'
    stdout, stderr = age_logs.spool_paths(archive, identifier)
    stdout = Path(str(stdout).replace('%j', identifier))
    stderr = Path(str(stderr).replace('%j', identifier))
    stdout.write_bytes(b'worker output\n' * 100)
    stderr.write_bytes(b'worker error\n')
    status = tmp_path / 'state' / f'{identifier}.status.json'
    write_json(status, {'job_id': identifier, 'controller_id': '99', 'rule': 'translate_cds',
                        'log_archive': str(archive), 'stdout': str(stdout), 'stderr': str(stderr)})
    return archive, status, stdout, stderr


def contents(archive):
    with tarfile.open(archive, 'r:gz') as packed:
        return {m.name: packed.extractfile(m).read() for m in packed}


def test_concurrent_successful_and_failed_workers_share_one_readable_archive(tmp_path):
    workers = [worker(tmp_path, str(n)) for n in range(101, 109)]
    failed = {'state': 'FAILED', 'source': 'qacct', 'exit_code': '152:1'}
    with ThreadPoolExecutor(max_workers=4) as pool:
        list(pool.map(lambda row: age_logs.archive_worker(row[1], failed if row == workers[-1] else COMPLETE), workers))
    archive = workers[0][0]
    age_logs.finalize(archive)
    size = archive.stat().st_size
    # Repeated status checks and finalization must not duplicate or rewrite logs.
    age_logs.archive_worker(workers[-1][1], failed)
    age_logs.finalize(archive)
    assert archive.stat().st_size == size
    result = contents(archive)
    assert len(result) == 24
    for _, status, stdout, stderr in workers:
        assert not stdout.exists() and not stderr.exists()
        identifier = json.loads(status.read_text())['job_id']
        assert result[f'{identifier}/stdout.txt'] == b'worker output\n' * 100
        assert result[f'{identifier}/stderr.txt'] == b'worker error\n'
    assert json.loads(result['108/job.json'])['outcome'] == failed
    assert len(list(archive.parent.iterdir())) == 3  # Archive, journal, lock; no per-worker logs.
    if shutil.which('tar'):
        check = subprocess.run(['tar', '-tzf', str(archive)], text=True, capture_output=True, check=True)
        assert len(check.stdout.splitlines()) == 24


def test_queued_or_unconfirmed_workers_are_never_archived(tmp_path):
    archive, status, stdout, stderr = worker(tmp_path)
    for outcome in ({'state': 'RUNNING', 'source': 'qstat'}, {'state': 'FAILED', 'source': 'qstat'},
                    {'state': 'COMPLETED', 'source': 'missing'}):
        age_logs.archive_worker(status, outcome)
    assert stdout.exists() and stderr.exists() and not archive.exists()


@pytest.mark.parametrize('failure_point', ['_member', '_verify'])
def test_failed_verification_retains_raw_logs_and_retries_partial_append(tmp_path, monkeypatch, failure_point):
    archive, status, stdout, stderr = worker(tmp_path)
    original = getattr(age_logs, failure_point)
    monkeypatch.setattr(age_logs, failure_point, lambda *args: (_ for _ in ()).throw(ValueError('injected verification failure')))
    with pytest.raises(ValueError, match='injected verification failure'):
        age_logs.archive_worker(status, COMPLETE)
    assert stdout.exists() and stderr.exists()
    assert not json.loads(status.read_text()).get('logs_archived')
    monkeypatch.setattr(age_logs, failure_point, original)
    age_logs.archive_worker(status, COMPLETE)
    age_logs.finalize(archive)
    assert len(contents(archive)) == 3
    assert not stdout.exists() and not stderr.exists()


def test_interrupted_retirement_recovers_without_duplicate_entries(tmp_path, monkeypatch):
    archive, status, stdout, stderr = worker(tmp_path)
    original = age_logs._retire
    monkeypatch.setattr(age_logs, '_retire', lambda *args: (_ for _ in ()).throw(OSError('injected interruption')))
    with pytest.raises(OSError, match='injected interruption'):
        age_logs.archive_worker(status, COMPLETE)
    assert stdout.exists() and stderr.exists()
    monkeypatch.setattr(age_logs, '_retire', original)
    age_logs.archive_worker(status, COMPLETE)
    age_logs.finalize(archive)
    assert len(contents(archive)) == 3
    assert not stdout.exists() and not stderr.exists()


def test_new_worker_can_append_after_previous_finalization(tmp_path):
    archive, first, _, _ = worker(tmp_path)
    age_logs.archive_worker(first, COMPLETE)
    age_logs.finalize(archive)
    _, second, _, _ = worker(tmp_path, '102')
    age_logs.archive_worker(second, COMPLETE)
    age_logs.finalize(archive)
    assert len(contents(archive)) == 6


def test_recovery_collects_only_confirmed_workers_in_the_selected_run(tmp_path, monkeypatch):
    archive, status, stdout, stderr = worker(tmp_path)
    other, _, other_stdout, _ = worker(tmp_path / 'other_run')
    import scheduler
    monkeypatch.setattr(scheduler, 'age_states', lambda ids: {identifier: COMPLETE for identifier in ids})
    age_logs.recover(archive.parent, tmp_path / 'state')
    assert not stdout.exists() and not stderr.exists()
    assert json.loads(archive.with_suffix('.json').read_text())['finalized']
    assert not other.exists() and other_stdout.exists()


def test_archive_metadata_distinguishes_missing_logs_from_empty_logs(tmp_path):
    archive, status, stdout, stderr = worker(tmp_path)
    stdout.write_bytes(b'')
    stderr.unlink()
    age_logs.archive_worker(status, COMPLETE)
    age_logs.finalize(archive)
    result = contents(archive)
    assert result['101/stdout.txt'] == b''
    assert json.loads(result['101/job.json'])['missing_streams'] == ['stderr']


@pytest.mark.parametrize('exit_code', [0, 7])
def test_controller_finalizes_logs_on_success_and_failure(tmp_path, exit_code):
    script = tmp_path / 'controller.py'
    script.write_text(f'''
import json, os, sys
from pathlib import Path
sys.path.insert(0, {str(ROOT / 'workflow/scripts')!r})
from common import write_json
from age_logs import archive_worker, spool_paths
archive = os.environ['PHENORADAR_AGE_LOG_ARCHIVE']
stdout, stderr = [Path(str(p).replace('%j', '101')) for p in spool_paths(archive, 'worker')]
stdout.write_text('complete output')
stderr.write_text('diagnostic')
status = Path({str(tmp_path / 'state/101.status.json')!r})
write_json(status, {{'job_id': '101', 'log_archive': archive, 'stdout': str(stdout), 'stderr': str(stderr)}})
archive_worker(status, {COMPLETE!r})
raise SystemExit({exit_code})
''')
    env = {**os.environ, 'JOB_ID': '99'}
    result = subprocess.run([sys.executable, str(ROOT / 'workflow/scripts/age_logs.py'), 'run',
                             '--usage', str(tmp_path / 'jobs/storage'), '--state', str(tmp_path / 'state'),
                             '--', sys.executable, str(script)], env=env, text=True, capture_output=True, timeout=30)
    assert result.returncode == exit_code, result.stdout + result.stderr
    archive = tmp_path / 'jobs/logs/age-workers-99.tar.gz'
    assert json.loads(archive.with_suffix('.json').read_text())['finalized']
    assert contents(archive)['101/stdout.txt'] == b'complete output'
