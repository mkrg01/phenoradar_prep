"""Public Slurm time settings preserve duration across execution backends."""
from pathlib import Path
import subprocess

import pytest
import yaml

from phase_config import array_concurrency, execution_settings, merge_slurm, resolve_array_size, time_minutes, validate_slurm, worker_budgets, workflow_jobs, write_profile

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize('duration,minutes', [
    ('1-00:00:00', 1440), ('3-00:00:00', 4320), ('14-00:00:00', 20160),
    ('02:00:00', 120), ('25:00:00', 1500), ('90', 90), ('90:30', 91),
    ('2-03', 3060), ('2-03:04', 3064), ('2-03:04:05', 3065),
    ('00:00:01', 1), ('00:01:01', 2),
])
def test_time_conversion_matches_slurm_minutes(duration, minutes):
    assert time_minutes(duration) == minutes


@pytest.mark.parametrize('duration', [
    None, True, 3600, 1.5, '', '-1', '0', '00:00:00', '1h',
    '1-00:00:00:00', '1--00', '01::00', 'infinite',
])
def test_invalid_or_unquoted_time_is_rejected(duration):
    with pytest.raises(ValueError, match='invalid Slurm time'):
        time_minutes(duration)


@pytest.mark.parametrize('build', [True, False])
def test_time_override_converts_without_changing_saved_resources(tmp_path, build):
    phase = 'build' if build else 'analysis'
    slurm = yaml.safe_load((ROOT / f'config/{phase}.yaml').read_text())['slurm']
    original = yaml.safe_dump(slurm)
    override = tmp_path / 'resources.yaml'
    override.write_text(yaml.safe_dump({'slurm': {
        'stages': {'controller': {'time': '2-00:00:00'}},
        'default_resources': {'time': '06:00:00'},
        'rules': {'example': {'time': '01:00:01'}},
    }}))
    settings = execution_settings(slurm, override, build=build)
    assert yaml.safe_dump(slurm) == original
    assert settings['stages']['controller']['time'] == '2-00:00:00'
    profile = write_profile(tmp_path / 'profile', settings)
    actual = yaml.safe_load((profile / 'config.yaml').read_text())
    assert actual['default-resources']['runtime'] == 360
    assert actual['default-resources']['mem_mb'] == 8000
    assert actual['set-resources']['example'] == {'runtime': 61}
    assert 'time' not in actual['default-resources']
    assert settings['rules']['example']['time'] == '01:00:01'


@pytest.mark.parametrize('section', ['stages', 'default_resources', 'rules'])
def test_all_time_fields_reject_yaml_numeric_values(section):
    slurm = validate_slurm(yaml.safe_load((ROOT / 'config/analysis.yaml').read_text())['slurm'], build=False)
    # YAML 1.1 can silently turn an unquoted 12:00:00 into the integer 43200.
    resource = yaml.safe_load('time: 12:00:00')
    if section == 'stages':
        slurm['stages']['controller'].update(resource)
    elif section == 'rules':
        slurm['rules']['example'] = resource
    else:
        slurm['default_resources'].update(resource)
    with pytest.raises(ValueError, match='invalid Slurm time'):
        validate_slurm(slurm, build=False)


@pytest.fixture
def public_slurm():
    return yaml.safe_load((ROOT / 'config/build.yaml').read_text())['slurm']


def test_build_config_has_no_array_size_setting(public_slurm):
    assert 'array_size' not in public_slurm
    assert 'array_size' not in validate_slurm(public_slurm)


def test_legacy_stage_allocations_become_one_sample_without_mutating_snapshot(public_slurm):
    import copy
    legacy = copy.deepcopy(public_slurm)
    legacy['per_job_resources'].pop('sample')
    legacy['per_job_resources'].update(
        assembly={'cpus': 8, 'mem_gb': 128, 'time': '14-00:00:00'},
        busco={'cpus': 4, 'mem_gb': 32, 'time': '1-00:00:00'},
        quant={'cpus': 16, 'mem_gb': 64, 'time': '02:00:00'})
    before = copy.deepcopy(legacy)
    normalized = validate_slurm(legacy)
    assert normalized['stages']['sample'] == {'cpus': 16, 'mem_gb': 128, 'time': '15-02:00:00'}
    assert set(normalized['stages']) == {'sample', 'controller'}
    assert legacy == before


def test_legacy_assembly_retry_sizes_whole_sample(public_slurm, tmp_path):
    override = tmp_path / 'retry.yaml'
    override.write_text('slurm:\n  per_job_resources:\n    assembly:\n      mem_gb: 256\n')
    assert execution_settings(public_slurm, override)['stages']['sample']['mem_gb'] == 256


@pytest.mark.parametrize('value', [None, 1, 1000])
def test_legacy_array_size_is_ignored_without_mutation(public_slurm, value):
    public_slurm['array_size'] = value
    assert 'array_size' not in validate_slurm(public_slurm)
    assert public_slurm['array_size'] == value


def test_finite_public_config_preserves_legacy_profile_and_array_caps(public_slurm, tmp_path):
    import copy
    public_slurm['total_limits']['jobs'] = 64
    legacy = validate_slurm(public_slurm)
    legacy.pop('total_limits')
    saved = copy.deepcopy(public_slurm)
    new_profile = write_profile(tmp_path / 'new', public_slurm)
    old_profile = write_profile(tmp_path / 'old', legacy)
    assert (new_profile / 'config.yaml').read_bytes() == (old_profile / 'config.yaml').read_bytes()
    assert public_slurm == saved
    assert legacy['default_resources'] == {'mem_gb': 8, 'time': '1-00:00:00'}
    for stage in ('assembly', 'busco', 'quant'):
        assert array_concurrency(legacy, 'sample') == 64


def test_new_override_updates_both_legacy_caps_without_mutation(public_slurm, tmp_path):
    import copy
    legacy = validate_slurm(public_slurm)
    legacy['concurrency'], legacy['jobs'] = 3, 7
    original = copy.deepcopy(legacy)
    retry = tmp_path / 'retry.yaml'
    retry.write_text(yaml.safe_dump({'slurm': {'total_limits': {'jobs': 5, 'mem_gb': 512},
                                             'per_job_resources': {'sample': {'cpus': 8}}}}))
    updated = execution_settings(legacy, retry)
    assert updated['concurrency'] == updated['jobs'] == 5
    assert updated['stages']['sample'] == {'cpus': 8, 'mem_gb': 128, 'time': '16-00:00:00'}
    assert legacy == original
    assert execution_settings(legacy)['concurrency'] == 3
    assert execution_settings(legacy)['jobs'] == 7


@pytest.mark.parametrize('limits,expected', [
    ({'jobs': 3}, 3),
    ({'cpus': 10}, 2),
    ({'mem_gb': 512}, 4),
    ({'cpus': 10, 'mem_gb': 128}, 1),
])
def test_array_limits_use_all_budgets(public_slurm, limits, expected):
    slurm = merge_slurm(public_slurm, {'total_limits': limits})
    assert array_concurrency(slurm, 'sample') == expected
    # Budgeting changes concurrency, never the job request.
    assert slurm['stages']['sample']['cpus'] == 4
    assert slurm['stages']['sample']['mem_gb'] == 128


@pytest.mark.parametrize('key,value', [('cpus', 3), ('mem_gb', 127)])
def test_single_array_job_must_fit_budget(public_slurm, key, value):
    slurm = merge_slurm(public_slurm, {'total_limits': {key: value}})
    with pytest.raises(ValueError, match='sample requires.*exceeding total_limits'):
        array_concurrency(slurm, 'sample')


def test_worker_profile_reserves_controller_and_preserves_requests(public_slurm, tmp_path):
    slurm = merge_slurm(public_slurm, {'total_limits': {'cpus': 33, 'mem_gb': 512}})
    profile = yaml.safe_load((write_profile(tmp_path / 'profile', slurm) / 'config.yaml').read_text())
    assert profile['resources'] == {'workflow_cpus': 32, 'workflow_mem_mb': 504000}
    assert profile['set-resource-scopes'] == {'workflow_cpus': 'global', 'workflow_mem_mb': 'global'}
    assert profile['set-threads']['odb_map'] == 16
    assert profile['set-resources']['odb_map']['mem_mb'] == 128000
    assert 'mem_mb' not in profile['resources']  # Do not let Snakemake clip the actual request.


@pytest.mark.parametrize('key,value', [('cpus', 1), ('mem_gb', 8)])
def test_controller_must_leave_room_for_workers(public_slurm, key, value):
    slurm = merge_slurm(public_slurm, {'total_limits': {key: value}})
    with pytest.raises(ValueError, match='leaves no worker budget'):
        worker_budgets(slurm)


@pytest.mark.parametrize('key,value', [('jobs', False), ('jobs', 0), ('cpus', 0),
                                      ('mem_gb', -1), ('cpus', True), ('mem_gb', 1.5)])
def test_invalid_total_limits(public_slurm, key, value):
    public_slurm['total_limits'][key] = value
    with pytest.raises(ValueError, match='must be (a positive integer or null|positive)'):
        validate_slurm(public_slurm)


def test_ambiguous_and_malformed_resources_rejected(public_slurm):
    public_slurm['jobs'] = 8
    with pytest.raises(ValueError, match='do not mix'):
        validate_slurm(public_slurm)
    del public_slurm['jobs']
    public_slurm['per_job_resources']['sample'] = None
    with pytest.raises(ValueError, match='invalid resources for sample'):
        validate_slurm(public_slurm)


@pytest.mark.parametrize('phase', ['build', 'analysis'])
def test_null_job_limit_is_explicitly_unlimited(phase, tmp_path):
    public = yaml.safe_load((ROOT / f'config/{phase}.yaml').read_text())['slurm']
    assert public['total_limits']['jobs'] is None
    slurm = validate_slurm(public, build=phase == 'build')
    assert slurm['jobs'] is None
    assert workflow_jobs(slurm) == 'unlimited'
    profile = yaml.safe_load((write_profile(tmp_path / 'profile', slurm) / 'config.yaml').read_text())
    assert profile['jobs'] == 'unlimited'
    if phase == 'build':
        assert array_concurrency(slurm, 'sample') is None


def test_null_override_removes_legacy_count_caps_but_keeps_cpu_budget(public_slurm, tmp_path):
    slurm = validate_slurm(public_slurm)
    slurm.update(jobs=2, concurrency=2)
    slurm['total_limits']['cpus'] = 16
    retry = tmp_path / 'retry.yaml'
    retry.write_text('slurm:\n  total_limits:\n    jobs: null\n')
    updated = execution_settings(slurm, retry)
    assert updated['jobs'] is updated['concurrency'] is None
    assert array_concurrency(updated, 'sample') == 4
    assert slurm['jobs'] == slurm['concurrency'] == 2


@pytest.mark.parametrize('configuration,expected', [
    ('MaxArraySize = 1001\nSchedulerParameters = (null)\n', 1000),
    ('MaxArraySize = 4000001\nSchedulerParameters = (null)\n', 4000000),
    ('MaxArraySize = 100001\nSchedulerParameters = bf_continue,max_array_tasks=1000\n', 1000),
    ('MaxArraySize = 3\nSchedulerParameters = max_array_tasks=10,bf_continue\n', 2),
    ('  MaxArraySize = 11\nSchedulerParameters = bf_continue, max_array_tasks=4\n', 4),
])
def test_detect_array_size_from_slurm(monkeypatch, configuration, expected):
    calls = []
    def query(command, **kwargs):
        calls.append(command)
        assert kwargs == {'text': True, 'stderr': subprocess.PIPE, 'timeout': 10}
        return configuration
    monkeypatch.setattr(subprocess, 'check_output', query)
    assert resolve_array_size() == expected
    assert calls == [['scontrol', 'show', 'config']]


@pytest.mark.parametrize('configuration', [
    '', 'SchedulerParameters = (null)\n', 'MaxArraySize = 1001\n',
    'MaxArraySize = invalid\nSchedulerParameters = (null)\n',
    'MaxArraySize = 0\nSchedulerParameters = (null)\n',
    'MaxArraySize = 1\nSchedulerParameters = (null)\n',
    'MaxArraySize = 1001\nSchedulerParameters = max_array_tasks=0\n',
    'MaxArraySize = 1001\nSchedulerParameters = max_array_tasks=invalid\n',
    'MaxArraySize = 1001\nSchedulerParameters = max_array_tasks\n',
])
def test_unusable_slurm_array_limits_fail_explicitly(monkeypatch, configuration):
    monkeypatch.setattr(subprocess, 'check_output', lambda *args, **kwargs: configuration)
    with pytest.raises(ValueError, match='Slurm'):
        resolve_array_size()


@pytest.mark.parametrize('error', [FileNotFoundError('scontrol'),
    subprocess.CalledProcessError(1, ['scontrol', 'show', 'config']),
    subprocess.TimeoutExpired(['scontrol', 'show', 'config'], 10)])
def test_array_limit_query_failure_is_reported(monkeypatch, error):
    def fail(*args, **kwargs):
        raise error
    monkeypatch.setattr(subprocess, 'check_output', fail)
    with pytest.raises(ValueError, match='Cannot query Slurm array limits'):
        resolve_array_size()
