"""Exercise the real Snakemake scheduler with synthetic jobs and global budgets."""
import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

from phase_config import merge_slurm, write_profile

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def budget_workflow(tmp_path):
    snakefile = tmp_path / 'Snakefile'
    snakefile.write_text('''
import sys
sys.path.insert(0, SCRIPTS)

rule all:
    input: expand("done/{n}", n=range(4))

rule work:
    output: "done/{n}"
    threads: 2
    resources: mem_mb=lambda wildcards: 2500
    shell: "PYTHON worker.py {output}"

from resource_limits import apply_workflow_limits
apply_workflow_limits(workflow)
'''.replace('SCRIPTS', repr(str(ROOT / 'workflow/scripts'))).replace('PYTHON', sys.executable))
    (tmp_path / 'worker.py').write_text('''
import fcntl, json, sys, time
from pathlib import Path

def event(action):
    with open('events.jsonl', 'a') as handle:
        fcntl.flock(handle, fcntl.LOCK_EX)
        handle.write(json.dumps([action, sys.argv[1]]) + '\\n')
        handle.flush()

event('start')
time.sleep(0.5)
Path(sys.argv[1]).parent.mkdir(exist_ok=True)
Path(sys.argv[1]).write_text('done')
event('end')
''')
    return tmp_path, snakefile


def run_workflow(project, overrides, *, dry_run=False):
    executable = shutil.which('snakemake')
    if not executable: pytest.skip('Snakemake required')
    path, snakefile = project
    public = yaml.safe_load((ROOT / 'config/analysis.yaml').read_text())['slurm']
    slurm = merge_slurm(public, overrides, build=False)
    profile = write_profile(path / 'profile', slurm)
    command = [executable, '--snakefile', str(snakefile), '--profile', str(profile)]
    # Dry-run the actual remote profile as well as running its resource accounting
    # through Snakemake's local scheduler, without submitting to a real cluster.
    if dry_run: command += ['--dry-run']
    else: command += ['--executor', 'local', '--cores', '8']
    return subprocess.run(command + ['--', 'all'], cwd=path, text=True, capture_output=True, timeout=90)


@pytest.mark.parametrize('limits,job,peak', [
    ({'cpus': 5}, {}, 2),
    ({'mem_gb': 13}, {}, 2),
    ({'cpus': 5}, {'cpus': 3}, 1),
    ({'mem_gb': 13}, {'mem_gb': 3}, 1),
])
def test_budget_bounds_real_parallel_execution(budget_workflow, limits, job, peak):
    result = run_workflow(budget_workflow, {'total_limits': limits, 'per_job_resources': {'work': job}})
    assert result.returncode == 0, result.stdout + result.stderr
    active = maximum = started = 0
    for action, _ in map(json.loads, (budget_workflow[0] / 'events.jsonl').read_text().splitlines()):
        active += 1 if action == 'start' else -1
        started += action == 'start'
        maximum = max(maximum, active)
    assert (active, started, maximum) == (0, 4, peak)


@pytest.mark.parametrize('limits,resource', [({'cpus': 2}, 'cpus'), ({'mem_gb': 10}, 'mem_gb')])
def test_scheduler_rejects_oversized_job_without_clipping(budget_workflow, limits, resource):
    result = run_workflow(budget_workflow, {'total_limits': limits})
    assert result.returncode != 0
    assert f'total_limits.{resource} budget' in result.stdout + result.stderr
    assert not (budget_workflow[0] / 'events.jsonl').exists()


def test_remote_profile_counts_dynamic_memory_and_thread_overrides(budget_workflow):
    result = run_workflow(budget_workflow, {'total_limits': {'cpus': 5, 'mem_gb': 13},
                                          'per_job_resources': {'work': {'cpus': 3}}}, dry_run=True)
    assert result.returncode == 0, result.stdout + result.stderr
    output = result.stdout + result.stderr
    assert 'workflow_cpus=3' in output
    assert 'workflow_mem_mb=<TBD>' in output  # Input-dependent accounting resolves before scheduling.
    assert 'mem_mb=2500' in output


def test_budget_checks_memory_after_input_is_generated(budget_workflow):
    path, snakefile = budget_workflow
    source = snakefile.read_text().replace('rule work:\n', 'rule work:\n    input: "memory.txt"\n')
    source = source.replace('mem_mb=lambda wildcards: 2500',
                            'mem_mb=lambda wildcards, input: int(open(input[0]).read())')
    source = source.replace('from resource_limits import apply_workflow_limits', '''
rule memory:
    output: "memory.txt"
    resources: mem_mb=1000
    shell: "echo 6000 > {output}"

from resource_limits import apply_workflow_limits''')
    snakefile.write_text(source)
    result = run_workflow(budget_workflow, {'total_limits': {'mem_gb': 13}})
    assert result.returncode != 0
    assert (path / 'memory.txt').read_text().strip() == '6000'
    assert 'work requires 6 mem_gb' in result.stdout + result.stderr
    assert not (path / 'events.jsonl').exists()
