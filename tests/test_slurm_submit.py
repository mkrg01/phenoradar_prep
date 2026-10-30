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


def test_gate_limits_distinct_actual_submissions_and_not_receipt_reuse(tmp_path):
    current = [100.0]
    starts = []
    def sleep(delay): current[0] += delay
    def run(argv, **kwargs):
        starts.append(current[0])
        return CompletedProcess(argv, 0, str(100 + len(starts)), '')
    for comment in ('one', 'two', 'three', 'three'):
        submit([*ARGS[:2], '--comment=' + comment, *ARGS[4:]], '/bin/sbatch', tmp_path,
               run, sleep, lambda: current[0])
    assert starts == [100.0, 102.0, 104.0]


def test_gate_serializes_rpc_even_if_first_submitter_pauses(tmp_path):
    import threading
    from concurrent.futures import ThreadPoolExecutor
    entered, release, second = threading.Event(), threading.Event(), threading.Event()
    current = [100.0]
    def run(argv, **kwargs):
        if not entered.is_set():
            entered.set()
            assert release.wait(3)
        else:
            second.set()
        return CompletedProcess(argv, 0, '123', '')
    def submit_one(comment):
        return submit([*ARGS[:2], '--comment=' + comment, *ARGS[4:]], '/bin/sbatch', tmp_path,
                      run, lambda delay: current.__setitem__(0, current[0] + delay), lambda: current[0])
    with ThreadPoolExecutor(max_workers=2) as pool:
        first = pool.submit(submit_one, 'one')
        assert entered.wait(3)
        next_job = pool.submit(submit_one, 'two')
        try:
            assert not second.wait(0.1)
        finally:
            release.set()
        assert first.result(3) == next_job.result(3) == '123'


@pytest.mark.parametrize('mode,executor,guarded,disabled', [
    ('default', 'slurm', True, True), ('remote', 'slurm', True, False),
    ('default', 'local', True, False), ('default', 'slurm', False, False)])
def test_only_guarded_controller_disables_scheduler_timers(monkeypatch, tmp_path, mode, executor, guarded, disabled):
    from types import SimpleNamespace
    from snakemake_interface_executor_plugins.settings import ExecMode
    import slurm_submit
    sentinel = object()
    workflow = SimpleNamespace(exec_mode=ExecMode.DEFAULT if mode == 'default' else ExecMode.REMOTE,
                               executor_plugin=SimpleNamespace(name=executor),
                               scheduler=SimpleNamespace(job_rate_limiter=sentinel))
    shim = Path(slurm_submit.__file__).resolve().parents[1] / 'slurm/sbatch'
    monkeypatch.setattr(slurm_submit.shutil, 'which', lambda _: str(shim if guarded else '/bin/sbatch'))
    monkeypatch.setenv('PHENORADAR_REAL_SBATCH', '/bin/sbatch')
    monkeypatch.setenv('PHENORADAR_SUBMIT_STATE', str(tmp_path))
    probes = []
    monkeypatch.setattr(slurm_submit.subprocess, 'run', lambda argv, **kw: probes.append(argv))
    slurm_submit.configure_scheduler(workflow)
    assert (workflow.scheduler.job_rate_limiter is None) == disabled
    assert len(probes) == int(disabled)


def test_unavailable_accounting_stops_before_disabling_limiter(monkeypatch, tmp_path):
    from types import SimpleNamespace
    from snakemake_interface_executor_plugins.settings import ExecMode
    from subprocess import CalledProcessError
    import slurm_submit
    sentinel = object()
    workflow = SimpleNamespace(exec_mode=ExecMode.DEFAULT, executor_plugin=SimpleNamespace(name='slurm'),
                               scheduler=SimpleNamespace(job_rate_limiter=sentinel))
    shim = Path(slurm_submit.__file__).resolve().parents[1] / 'slurm/sbatch'
    monkeypatch.setattr(slurm_submit.shutil, 'which', lambda _: str(shim))
    monkeypatch.setenv('PHENORADAR_REAL_SBATCH', '/bin/sbatch')
    monkeypatch.setenv('PHENORADAR_SUBMIT_STATE', str(tmp_path))
    def fail(*args, **kw): raise CalledProcessError(1, 'sacct')
    monkeypatch.setattr(slurm_submit.subprocess, 'run', fail)
    with pytest.raises(RuntimeError, match='accounting.*required'):
        slurm_submit.configure_scheduler(workflow)
    assert workflow.scheduler.job_rate_limiter is sentinel


def test_exhausted_worker_slots_wait_without_reevaluation_timers(monkeypatch, tmp_path):
    """Exercise the pinned scheduler after completed workers have filled its slots."""
    import sys
    import threading
    from types import SimpleNamespace
    import snakemake.api
    from snakemake.scheduler import JobScheduler, JobRateLimiter
    from snakemake.settings.types import MaxJobsPerTimespan
    from snakemake_interface_executor_plugins.settings import ExecMode
    import slurm_submit

    stopped, waiting = threading.Event(), threading.Event()
    class EndProbe(Exception): pass
    class Wakeups:
        def __init__(self): self.semaphore, self.calls = threading.Semaphore(1), 0
        def acquire(self):
            self.calls += 1
            if self.calls == 2: waiting.set()
            self.semaphore.acquire()
        def release(self): self.semaphore.release()
    class Job:
        priority = 0
        scheduler_resources = {'_cores': 1, '_nodes': 1}
        def reset_params_and_resources(self): pass
        async def inputsize(self): return 1
    class Dag:
        queue_input_jobs = []
        @property
        def ready_jobs(self):
            if stopped.is_set(): raise EndProbe
            return {job}
        def register_running(self, jobs): assert not jobs
        async def temp_size(self, job): return 0
    job = Job()
    scheduler = object.__new__(JobScheduler)
    workflow = SimpleNamespace(exec_mode=ExecMode.DEFAULT, executor_plugin=SimpleNamespace(name='slurm'),
                               scheduler=scheduler, dag=Dag(),
                               remote_execution_settings=SimpleNamespace(immediate_submit=False))
    scheduler.workflow = workflow
    scheduler.dryrun = scheduler.touch = scheduler.keepgoing = False
    scheduler.running, scheduler.failed = {object()}, set()
    scheduler._tofinish, scheduler._toerror = [], []
    scheduler._errors = False
    scheduler._executor_error = scheduler._user_kill = scheduler._run_performed = None
    scheduler._last_job_selection_empty, scheduler.greediness, scheduler.subsample = True, 1.0, None
    scheduler._executor = SimpleNamespace(cancel=lambda: None)
    scheduler.global_resources = {'_cores': 1, '_nodes': 1, '_job_count': sys.maxsize}
    scheduler.resources = dict(scheduler.global_resources, _nodes=0)
    scheduler._job_selector = scheduler.job_selector_greedy
    scheduler._lock, scheduler._open_jobs = threading.Lock(), Wakeups()
    scheduler.job_rate_limiter = JobRateLimiter(MaxJobsPerTimespan(1, '2s'))
    timers = []
    scheduler._schedule_reevalutation = lambda delay: timers.append(delay)
    shim = Path(slurm_submit.__file__).resolve().parents[1] / 'slurm/sbatch'
    monkeypatch.setattr(slurm_submit.shutil, 'which', lambda _: str(shim))
    monkeypatch.setenv('PHENORADAR_REAL_SBATCH', '/bin/sbatch')
    monkeypatch.setenv('PHENORADAR_SUBMIT_STATE', str(tmp_path))
    monkeypatch.setattr(slurm_submit.subprocess, 'run', lambda *args, **kw: None)
    slurm_submit.configure_scheduler(workflow)
    errors = []
    def probe():
        try: scheduler.schedule()
        except EndProbe: pass
        except Exception as error: errors.append(error)
    thread = threading.Thread(target=probe, daemon=True)
    thread.start()
    try:
        assert waiting.wait(3)
        assert timers == []
        assert scheduler._open_jobs.calls == 2
    finally:
        stopped.set(); scheduler._open_jobs.release(); thread.join(3)
    assert not thread.is_alive() and errors == []
