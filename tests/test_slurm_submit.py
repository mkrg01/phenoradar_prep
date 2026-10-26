import json
from pathlib import Path
from subprocess import CompletedProcess

import pytest

from slurm_submit import submit, UnknownSubmission, check_analysis_submissions


ARGS = ['--parsable', '--job-name=workflow-uuid', '--comment', 'rule_batch_001', '--chdir=/tmp']


def test_timeout_recovers_accepted_job_without_resubmission(tmp_path):
    calls = []
    token = None
    def run(argv, **kwargs):
        nonlocal token
        calls.append(argv)
        if argv[0] == '/bin/sbatch':
            token = next(a.split('phenoradar:')[1] for a in argv if a.startswith('--comment='))
            return CompletedProcess(argv, 1, '', 'Socket timed out on send/recv operation')
        return CompletedProcess(argv, 0, f'123|rule_batch_001;phenoradar:{token}|/tmp\n', '')
    assert submit(ARGS, '/bin/sbatch', tmp_path, run, lambda _: None) == '123'
    assert sum(c[0] == '/bin/sbatch' for c in calls) == 1
    assert submit(ARGS, '/bin/sbatch', tmp_path, run, lambda _: None) == '123'
    assert sum(c[0] == '/bin/sbatch' for c in calls) == 1


def test_unknown_timeout_never_blindly_resubmits_even_on_retry(tmp_path):
    calls = []
    def run(argv, **kwargs):
        calls.append(argv)
        return CompletedProcess(argv, 1 if argv[0] == '/bin/sbatch' else 0,
                                '', 'Socket timed out on send/recv operation')
    for argv in (ARGS, ARGS + ['--wrap=attempt2']):
        with pytest.raises(UnknownSubmission):
            submit(argv, '/bin/sbatch', tmp_path, run, lambda _: None)
    assert sum(c[0] == '/bin/sbatch' for c in calls) == 1


def test_confirmed_temporary_rejection_retries_and_preserves_script_position(tmp_path):
    calls = []
    def run(argv, **kwargs):
        calls.append(argv)
        return CompletedProcess(argv, int(len(calls) < 3), '456\n' if len(calls) == 3 else '',
                                'Slurm temporarily unable to accept job')
    assert submit(ARGS + ['script.sh'], '/bin/sbatch', tmp_path, run, lambda _: None) == '456'
    assert len(calls) == 3
    assert all(c[-1] == 'script.sh' and c[1].startswith('--comment=') for c in calls)


def test_restart_blocks_only_surviving_or_ambiguous_workers_of_this_analysis(tmp_path):
    analysis = tmp_path / 'analysis'
    receipt = tmp_path / 'worker.json'
    receipt.write_text(json.dumps({'analysis': str(analysis), 'state': 'accepted', 'job_id': '123'}))
    with pytest.raises(ValueError, match='queued/running worker'):
        check_analysis_submissions(tmp_path, analysis, {'123'})
    check_analysis_submissions(tmp_path, analysis, set())
    receipt.write_text(json.dumps({'analysis': str(analysis), 'state': 'unknown'}))
    with pytest.raises(ValueError, match='unresolved worker'):
        check_analysis_submissions(tmp_path, analysis, set())
    check_analysis_submissions(tmp_path, tmp_path / 'other-analysis', set())
