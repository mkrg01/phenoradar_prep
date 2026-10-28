"""Execution observations must stay cheap and distinguish evidence from absence."""
import json
import logging
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

import analysis
from analysis_runtime import ProgressHandler, RunState, read_state
from analysis_status import observe_slurm, required_outputs, status
from common import write_json


@pytest.fixture
def saved(tmp_path):
    path = tmp_path / 'results/test'
    path.mkdir(parents=True)
    write_json(path / 'analysis.json', {'name': 'test', 'root': str(tmp_path),
               'build': {'path': 'database/manifest.json'}, 'selection': [], 'pipeline': {}})
    return path


def receipt(path, **fields):
    saved = dict(created_at='2026-10-26T19:50:00+00:00', target='all',
                 state='submitted', job_id='123', **fields)
    write_json(path / 'jobs/submission_0003.json', saved)
    return saved


def queue(root, state=None, reason='None', accounting=None):
    def run(command, **kwargs):
        if command[0] == 'squeue':
            text = f'123|{state}|{reason}|{root}\n' if state else ''
        else:
            text = f'123|{accounting}|0:0|{root}|\n' if accounting else ''
        return subprocess.CompletedProcess(command, 0, text, '')
    return run


def complete(path, when='2026-10-27T00:00:00+00:00'):
    for name in required_outputs({}):
        target = path / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text('data\n')
    write_json(path / 'completed.json', {'created_at': when, 'files': []})


def test_default_status_never_verifies_or_collects_storage(saved, monkeypatch):
    def forbidden(*args, **kwargs):
        pytest.fail('heavy verification was called by default status')
    monkeypatch.setattr(analysis, 'load', forbidden)
    monkeypatch.setattr(analysis, 'verify', forbidden)
    before = {p: p.stat().st_mtime_ns for p in saved.rglob('*')}
    report = analysis.status(saved)
    assert report['state'] == 'prepared'
    assert report['verification']['performed'] is False
    assert report['storage']['state'] == 'not_checked'
    assert {p: p.stat().st_mtime_ns for p in saved.rglob('*')} == before


@pytest.mark.parametrize('state,expected', [('PENDING', 'pending'), ('RUNNING', 'running'),
    ('OUT_OF_MEMORY', 'failed'), ('TIMEOUT', 'failed'), ('CANCELLED', 'failed')])
def test_slurm_state_and_reason(saved, state, expected):
    receipt(saved)
    report = status(saved, run=queue(saved.parent.parent, state, 'Resources'))
    assert report['state'] == expected
    assert report['controller']['reason'] == 'Resources'


def test_queue_disappearance_and_unavailable_accounting_do_not_mean_success(saved):
    receipt(saved)
    assert status(saved, run=queue(saved.parent.parent))['state'] == 'unknown'
    def unavailable(*args, **kwargs):
        raise FileNotFoundError('squeue/sacct unavailable')
    report = status(saved, run=unavailable)
    assert report['state'] == 'unknown'
    assert len(report['controller']['errors']) == 2


def test_impossible_dependency_requires_action(saved):
    receipt(saved)
    report = status(saved, run=queue(saved.parent.parent, 'PENDING', 'DependencyNeverSatisfied'))
    assert report['state'] == 'blocked'


def test_completed_controller_needs_current_record_and_outputs(saved):
    receipt(saved)
    observed = queue(saved.parent.parent, accounting='COMPLETED')
    assert status(saved, run=observed)['state'] == 'incomplete'
    complete(saved, when='2026-10-25T00:00:00+00:00')
    assert status(saved, run=observed)['state'] == 'incomplete'
    complete(saved)
    assert status(saved, run=observed)['state'] == 'complete'
    (saved / 'phenoradar_inputs/tpm.tsv').unlink()
    report = status(saved, run=observed)
    assert report['state'] == 'incomplete'
    assert 'phenoradar_inputs/tpm.tsv' in report['completion']['missing']


def test_runtime_phase_100_percent_does_not_mean_whole_completion(saved, monkeypatch):
    receipt(saved)
    monkeypatch.setenv('SLURM_JOB_ID', '123')
    state = RunState(saved, 'all', saved / 'jobs/submission_0003.resources.json')
    state.phase('exporting', progress={'done': 824, 'total': 824})
    report = status(saved, run=queue(saved.parent.parent, 'RUNNING'))
    assert report['state'] == 'finalizing'
    assert report['phase'] == 'exporting'
    assert not report['completion']['current']
    state.phase('failed', last_error='export failed')
    assert status(saved, run=queue(saved.parent.parent))['state'] == 'failed'


def test_stale_runtime_and_resource_json_are_not_current_runs(saved):
    receipt(saved)
    write_json(saved / 'jobs/submission_0003.resources.json', {'slurm': {}})
    write_json(saved / 'jobs/run_state.json', {'controller_id': '122', 'phase': 'complete'})
    report = status(saved, run=queue(saved.parent.parent))
    assert report['runtime'] == {} and len(report['submissions']) == 1
    assert report['state'] == 'unknown'


def test_partial_target_is_not_whole_analysis_completion(saved):
    item = receipt(saved)
    item['target'] = 'alignments'
    write_json(saved / 'jobs/submission_0003.json', item)
    complete(saved)
    report = status(saved, run=queue(saved.parent.parent, accounting='COMPLETED'))
    assert report['state'] == 'target_complete'
    assert not report['completion']['required_outputs_present']


def test_slurm_only_accepts_known_ids_with_repository_workdir(tmp_path):
    commands = []
    def run(command, **kwargs):
        commands.append(command)
        if command[0] == 'squeue':
            text = f'123|RUNNING|None|{tmp_path}/other\n124|RUNNING|None|{tmp_path}\n125|RUNNING|None|{tmp_path}\n'
        else:
            text = f'123|FAILED|1:0|{tmp_path}/other|\n123.batch|FAILED|1:0|{tmp_path}|\n'
        return subprocess.CompletedProcess(command, 0, text, '')
    report = observe_slurm('123', tmp_path, '2026-10-26T00:00:00Z', ['124'], run)
    assert report['state'] == 'unknown'
    assert [w['job_id'] for w in report['workers']] == ['124']
    assert commands[1][commands[1].index('--jobs') + 1] == '123'
    assert '--starttime' in commands[1]


def test_explicit_verification_checks_recorded_outputs(saved, monkeypatch):
    complete(saved)
    artifact = saved / 'run.json'
    from dataset_assets import record
    write_json(saved / 'completed.json', {'created_at': '2026-10-27T00:00:00Z', 'files': [record(artifact)]})
    monkeypatch.setattr(analysis, 'load', lambda path: json.loads((saved / 'analysis.json').read_text()))
    monkeypatch.setattr('storage_management.storage_report', lambda *args: {'state': 'retained'})
    assert analysis.status(saved, verify_outputs=True)['verification']['performed']
    artifact.write_text('corrupt')
    with pytest.raises(ValueError, match='changed'):
        analysis.status(saved, verify_outputs=True)


def test_status_cli_json_is_read_only_even_without_snakemake(saved):
    cli = saved.parent.parent / 'cli'
    scripts = cli / 'workflow/scripts'
    shutil.copytree(Path(analysis.__file__).parent, scripts, ignore=shutil.ignore_patterns('__pycache__'))
    # Run a separate process: the normal entry point must not open SQLite.
    result = subprocess.run([sys.executable, str(scripts / 'analysis.py'), 'status', '--analysis', str(saved), '--json'],
                            capture_output=True, text=True, check=True)
    assert json.loads(result.stdout)['state'] == 'prepared'
    assert not (saved / 'jobs').exists()
    assert not (cli / '.cache').exists()
    assert not list(cli.rglob('*.pyc'))


def test_progress_final_write_and_worker_guard(saved, monkeypatch):
    monkeypatch.setenv('SLURM_JOB_ID', '123')
    state = RunState(saved, 'all')
    state.phase('computing')
    handler = ProgressHandler(state.path, state.token)
    def event(done, total):
        record = logging.LogRecord('snakemake', logging.INFO, '', 0, '', (), None)
        record.event, record.done, record.total = 'progress', done, total
        handler.emit(record)
    event(1, 10)
    event(2, 10)
    assert read_state(state.path)['progress']['done'] == 1
    event(10, 10)
    assert read_state(state.path)['progress']['done'] == 10
    handler.emit(logging.LogRecord('snakemake', logging.INFO, '', 0,
                 'Job 1 has been submitted with SLURM jobid 124', (), None))
    assert read_state(state.path)['recent_worker_ids'] == ['124']
    monkeypatch.setenv('SLURM_JOB_ID', '124')
    event(20, 20)
    assert read_state(state.path)['progress']['done'] == 10
    state.phase('exporting')
    assert 'progress' not in read_state(state.path)
