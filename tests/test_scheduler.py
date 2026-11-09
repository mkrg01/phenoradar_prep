"""Portable submission uses allocation evidence and verifies AGE accounting."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest
import yaml

import scheduler
from phase_config import validate_slurm, write_profile

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def clean_detection(monkeypatch):
    for key in ('PHENORADAR_SCHEDULER', 'SLURM_JOB_ID', 'JOB_ID', 'SGE_ROOT', 'NSLOTS'):
        monkeypatch.delenv(key, raising=False)
    scheduler.age_capabilities.cache_clear()
    yield
    scheduler.age_capabilities.cache_clear()


@pytest.mark.parametrize('commands,expected', [
    ({'sbatch', 'squeue', 'sacct'}, 'slurm'), ({'qsub', 'qstat', 'qacct'}, 'age'), (set(), None),
])
def test_detection_from_complete_client_set(clean_detection, monkeypatch, commands, expected):
    monkeypatch.setattr(shutil, 'which', lambda command: '/bin/' + command if command in commands else None)
    assert scheduler.detect(required=False) == expected


@pytest.mark.parametrize('backend,environment', [('slurm', {'SLURM_JOB_ID': '12'}),
                                               ('age', {'JOB_ID': '12', 'NSLOTS': '4'})])
def test_allocation_resolves_hosts_with_both_client_sets(clean_detection, monkeypatch, backend, environment):
    monkeypatch.setattr(shutil, 'which', lambda name: '/bin/' + name)
    monkeypatch.setenv(next(iter(environment)), next(iter(environment.values())))
    for key, value in environment.items():
        monkeypatch.setenv(key, value)
    assert scheduler.detect() == backend


def test_ambiguous_or_missing_scheduler_requires_explicit_override(clean_detection, monkeypatch):
    monkeypatch.setattr(shutil, 'which', lambda name: '/bin/' + name)
    with pytest.raises(ValueError, match='Both Slurm and AGE'):
        scheduler.detect()
    monkeypatch.setenv('PHENORADAR_SCHEDULER', 'age')
    assert scheduler.detect() == 'age'
    monkeypatch.delenv('PHENORADAR_SCHEDULER')
    monkeypatch.setattr(shutil, 'which', lambda name: None)
    with pytest.raises(ValueError, match='No scheduler'):
        scheduler.detect()


def test_age_discovers_single_node_pe_and_memory_resource(clean_detection, monkeypatch):
    monkeypatch.setenv('PHENORADAR_SCHEDULER', 'age')
    replies = {'-spl': 'mpi\ndef_slot\n', '-sc': '# name shortcut type\ns_vmem sv MEMORY\n',
               '-sp': 'allocation_rule $pe_slots\n'}
    monkeypatch.setattr(subprocess, 'check_output', lambda command, **kwargs: replies[command[1]])
    actual = scheduler.resolve({'scheduler': 'auto', 'age': {'parallel_environment': None}})
    assert actual['age']['parallel_environment'] == 'def_slot'
    assert actual['age']['memory_resource'] == 's_vmem'
    scheduler.age_capabilities.cache_clear()
    replies['-sp'] = 'allocation_rule $round_robin\n'
    with pytest.raises(ValueError, match='one node'):
        scheduler.resolve({'scheduler': 'auto'})


def age_execution():
    return {'scheduler': 'age', 'age': {'parallel_environment': 'def_slot', 'memory_resource': 's_vmem'}}


def test_shared_command_builder_preserves_units_and_dependency_semantics(tmp_path):
    job = {'cpus': 4, 'mem_gb': 128, 'time': '16-00:00:00'}
    kwargs = dict(script=tmp_path / 'job.sh', indices=[1, 3], concurrency=2, afterok=['10'], afterany=['11'])
    cmd = scheduler.batch_command(age_execution(), job, tmp_path, 'sample', 'logs/%A_%a.out', 'logs/%A_%a.err', **kwargs)
    assert 's_vmem=32000M,h_rt=1382400' in cmd
    assert cmd[cmd.index('-t') + 1] == '1-2' and cmd[cmd.index('-tc') + 1] == '2'
    assert cmd[cmd.index('-hold_jid') + 1] == '10,11'
    assert '$JOB_ID_$TASK_ID.out' in cmd[cmd.index('-o') + 1]
    slurm = scheduler.batch_command({'scheduler': 'slurm'}, job, tmp_path, 'sample', 'out', 'err', **kwargs)
    assert '--mem=128000M' in slurm and '--array=1,3%2' in slurm
    assert '--dependency=afterany:11,afterok:10' in slurm


def test_accounting_distinguishes_live_failed_and_unknown_tasks():
    commands = []
    def query(command, **kwargs):
        commands.append(command)
        if command[0] == 'qstat':
            output = '<job_info><job_list><JB_job_number>10</JB_job_number><state>r</state><tasks>1-2:1</tasks></job_list></job_info>'
            return subprocess.CompletedProcess(command, 0, output, '')
        if command[-1] == '10':
            return subprocess.CompletedProcess(command, 0, '====\njobnumber 10\ntaskid 3\nfailed 0\nexit_status 9\n', '')
        return subprocess.CompletedProcess(command, 1, '', 'accounting pending')
    actual = scheduler.age_states(['10_1', '10_2', '10_3', '20'], run=query)
    assert actual['10_1']['state'] == actual['10_2']['state'] == 'RUNNING'
    assert actual['10_3']['state'] == 'FAILED' and '20' not in actual
    assert commands.count(['qacct', '-j', '10']) == 1


def test_age_allocation_does_not_multiply_total_memory_twice(monkeypatch):
    monkeypatch.setenv('NSLOTS', '4')
    monkeypatch.setenv('PHENORADAR_TASK_MEM_MB', '128000')
    assert scheduler.allocation() == (4, 130217)


def test_age_profile_reuses_rule_resources_and_common_storage(tmp_path):
    public = yaml.safe_load((ROOT / 'config/build.yaml').read_text())['execution']
    settings = validate_slurm(public)
    settings.update(age_execution())
    profile = write_profile(tmp_path / 'profile', settings, root=tmp_path,
                            pipeline=tmp_path / 'pipeline.yaml', usage=tmp_path / 'jobs/storage')
    actual = yaml.safe_load((profile / 'config.yaml').read_text())
    assert actual['executor'] == 'cluster-generic'
    assert actual['set-resources']['odb_map']['mem_mb'] == 128000
    assert 'status-job --state' in actual['cluster-generic-status-cmd']
    assert not any(key.startswith('slurm-') for key in actual)


def test_slurm_storage_wrapper_runs_after_node_allocation(monkeypatch):
    monkeypatch.setenv('PHENORADAR_STORAGE_CONFIG', '/project/pipeline.yaml')
    monkeypatch.setenv('PHENORADAR_PROJECT_ROOT', '/project')
    monkeypatch.setenv('PHENORADAR_USAGE_ROOT', '/project/jobs/storage')
    actual = scheduler.wrap_slurm(['--job-name=test', '--wrap=python worker.py'])
    assert 'storage_runtime.py run' in actual[-1]
    assert 'bash -c' in actual[-1] and 'python worker.py' in actual[-1]


def test_age_unknown_worker_receipt_blocks_only_its_owner(tmp_path):
    from common import write_json
    directory = tmp_path / '.snakemake/cluster_submissions'
    write_json(directory / 'unknown.json', {'state': 'unknown', 'analysis': '/other'})
    scheduler.worker_guard(tmp_path, '/mine')
    write_json(directory / 'mine.json', {'state': 'unknown', 'analysis': '/mine'})
    with pytest.raises(ValueError, match='Unresolved worker'):
        scheduler.worker_guard(tmp_path, '/mine')


def test_accounting_delay_has_bounded_grace_and_never_implies_success(tmp_path, monkeypatch, capsys):
    from common import write_json
    path = tmp_path / '12.status.json'
    write_json(path, {'job_id': '12', 'missing_since': None})
    monkeypatch.setattr(sys, 'argv', ['scheduler.py', 'status-job', '--state', str(tmp_path), '12'])
    monkeypatch.setattr(scheduler, 'age_states', lambda ids: {})
    monkeypatch.setattr(scheduler.time, 'time', lambda: 100)
    scheduler.main()
    assert capsys.readouterr().out.strip() == 'running'
    monkeypatch.setattr(scheduler.time, 'time', lambda: 701)
    with pytest.raises(ValueError, match='completion is unconfirmed'):
        scheduler.main()
    monkeypatch.setattr(scheduler, 'age_states', lambda ids: {'12': {'state': 'COMPLETED'}})
    scheduler.main()
    assert capsys.readouterr().out.strip() == 'success'
    assert json.loads(path.read_text())['missing_since'] is None


def test_slurm_stdin_submission_uses_the_same_storage_wrapper(tmp_path, monkeypatch):
    from slurm_submit import submit
    monkeypatch.setenv('PHENORADAR_STORAGE_CONFIG', '/project/pipeline.yaml')
    monkeypatch.setenv('PHENORADAR_PROJECT_ROOT', '/project')
    monkeypatch.setenv('PHENORADAR_USAGE_ROOT', '/project/jobs/storage')
    calls = []
    def run(command, **kwargs):
        calls.append((command, kwargs))
        return subprocess.CompletedProcess(command, 0, '123\n', '')
    assert submit(['--parsable', '-'], '/bin/sbatch', tmp_path, run=run, stdin='echo result') == '123'
    script = calls[-1][1]['input']
    assert 'storage_runtime.py run' in script and 'bash -c' in script and 'echo result' in script


@pytest.mark.parametrize('allocated', [False, True])
def test_existing_pipeline_launcher_auto_detects_age_and_runs_a_real_dag(tmp_path, allocated):
    """Real Snakemake and common storage; only AGE's external service is doubled."""
    snakemake = Path(sys.executable).parent / 'snakemake'
    if not snakemake.exists():
        pytest.skip('Snakemake required')
    root = tmp_path / 'project'
    shutil.copytree(ROOT / 'workflow/scripts', root / 'workflow/scripts', ignore=shutil.ignore_patterns('__pycache__'))
    shutil.copy2(ROOT / 'run_pipeline.sh', root / 'run_pipeline.sh')
    (root / 'config').mkdir()
    (root / 'config/analysis.yaml').write_text(yaml.safe_dump({'execution': {
        'scheduler': 'auto', 'partition': None, 'total_limits': {'jobs': 2, 'cpus': 4, 'mem_gb': 16},
        'per_job_resources': {'controller': {'cpus': 1, 'mem_gb': 1, 'time': '01:00:00'}}}}))
    (root / 'pipeline.yaml').write_text('storage:\n  scratch_root: workspace\n')
    (root / 'workflow/Snakefile').write_text('''
rule all:
    input: 'result.txt'
rule first:
    output: 'first.txt'
    threads: 2
    resources: mem_mb=3000, runtime=3
    shell: "printf 'portable' > {output}"
rule second:
    input: 'first.txt'
    output: 'result.txt'
    resources: mem_mb=1000, runtime=1
    shell: 'cat {input} > {output}'
''')
    binary = tmp_path / 'age-bin'
    binary.mkdir()
    source = f'#!{sys.executable}\n' + '''import json,os,subprocess,sys
from pathlib import Path
name=Path(sys.argv[0]).name
args=sys.argv[1:]
state=Path(os.environ['AGE_TEST_STATE'])
state.mkdir(exist_ok=True)
if name=='qconf':
    print({'-spl':'def_slot','-sc':'s_vmem sv MEMORY','-sp':'allocation_rule $pe_slots','-sconf':'max_aj_tasks 1000'}[args[0]])
elif name=='qstat':
    print('<job_info/>')
elif name=='qacct':
    ident=args[args.index('-j')+1]
    entry=state/(ident+'.json')
    if not entry.exists(): sys.exit(1)
    code=json.loads(entry.read_text())['exit_status']
    print('====\\njobnumber '+ident+'\\ntaskid undefined\\nfailed 0\\nexit_status '+str(code))
elif name=='qsub':
    ident=str(100+len(list(state.glob('*.json'))))
    env=dict(os.environ, JOB_ID=ident, SGE_ROOT='/test/gridengine', NSLOTS=args[args.index('-pe')+2])
    if '-v' in args:
        env.update(field.split('=',1) for field in args[args.index('-v')+1].split(','))
    def log(flag): return args[args.index(flag)+1].replace('$JOB_ID',ident).replace('$TASK_ID','undefined')
    with open(log('-o'),'w') as out, open(log('-e'),'w') as err:
        result=subprocess.run(['bash',args[-1]], env=env, stdout=out, stderr=err)
    (state/(ident+'.json')).write_text(json.dumps({'exit_status':result.returncode,'arguments':args}))
    print(ident)
'''
    for name in ('qconf', 'qsub', 'qstat', 'qacct'):
        executable = binary / name
        executable.write_text(source)
        executable.chmod(0o755)
    env = {key: value for key, value in os.environ.items() if not key.startswith(('PHENORADAR_', 'SLURM_', 'SGE_'))
           and key not in {'JOB_ID', 'NSLOTS'}}
    env.update(PATH=str(binary) + os.pathsep + os.environ['PATH'],
               WORKFLOW_PYTHON=sys.executable, SNAKEMAKE_BIN=str(snakemake), AGE_TEST_STATE=str(tmp_path / 'state'))
    if allocated:
        env.update(JOB_ID='99', NSLOTS='4', SGE_ROOT='/test/gridengine', PHENORADAR_TASK_MEM_MB='16000')
    result = subprocess.run(['bash', str(root / 'run_pipeline.sh'), '--configfile', 'pipeline.yaml',
                             '--seconds-between-status-checks', '1', '--', 'all'],
                            cwd=root, env=env, capture_output=True, text=True, timeout=120)
    assert result.returncode == 0, result.stdout + result.stderr
    assert (root / 'result.txt').read_text() == 'portable'
    submissions = [json.loads(path.read_text()) for path in (tmp_path / 'state').glob('*.json')]
    if allocated:
        assert not submissions and 'localrule first' in result.stderr
    else:
        assert len(submissions) == 2 and all(entry['exit_status'] == 0 for entry in submissions)
        assert any('s_vmem=1500M,h_rt=180' in entry['arguments'] for entry in submissions)
    assert len(list((root / '.snakemake/storage').glob('*.json'))) == (1 if allocated else 2)
    assert not list((root / 'work/scratch').iterdir())


@pytest.mark.parametrize('failed,exit_status,expected', [
    ('0', '0', 'COMPLETED'),
    ('100 : assumedly after job', '152', 'FAILED'),
    ('0', '1', 'FAILED'),
])
def test_age_accounting_accepts_padded_field_values(failed, exit_status, expected):
    def query(command, **kwargs):
        if command[0] == 'qstat':
            output = '<job_info><queue_info/><job_info/></job_info>'
        else:
            output = ('==============================================================\n'
                      'jobnumber                10           \n'
                      'taskid                   undefined\n'
                      f'failed                   {failed}    \n'
                      f'exit_status              {exit_status}                  \n')
        return subprocess.CompletedProcess(command, 0, output, '')
    result = scheduler.age_states(['10'], run=query)
    assert result['10']['state'] == expected
    assert result['10']['source'] == 'qacct'
