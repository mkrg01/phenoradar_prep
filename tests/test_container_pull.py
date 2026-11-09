"""Compression limits, atomic publication and offline image reuse."""
import json
from pathlib import Path
import subprocess

import pytest

import container_pull as module


@pytest.mark.parametrize('env,cpus,memory', [({}, 8, 3814),
    ({'SLURM_JOB_ID': '12', 'SLURM_CPUS_PER_TASK': '2', 'SLURM_MEM_PER_NODE': '8192'}, 2, 2048),
    ({'SLURM_JOB_ID': '12', 'SLURM_CPUS_PER_TASK': '2', 'SLURM_MEM_PER_CPU': '4096'}, 2, 2048)])
def test_compression_observes_cpu_and_memory_allocation(env, cpus, memory):
    assert module.budgets(8, 16000, env) == (cpus, memory)


@pytest.mark.parametrize('cpus,memory,env', [(0, 16000, {}), (True, 16000, {}),
    (1, 1000, {}), (1, 16000, {'SLURM_JOB_ID': '12', 'SLURM_MEM_PER_NODE': '0'})])
def test_invalid_or_unlimited_budgets_fail_before_pull(cpus, memory, env):
    with pytest.raises(ValueError):
        module.budgets(cpus, memory, env)


@pytest.fixture
def runtime(tmp_path, monkeypatch):
    config = tmp_path / 'site.conf'
    config.write_text('allow setuid = no\nmksquashfs procs = 0\n# mksquashfs mem = 1G\n')
    state = {'suid': '0', 'fail': False, 'inspections': [], 'pulls': [], 'environments': []}
    def output(command, **kwargs):
        if command[1:] == ['buildcfg']:
            return f'APPTAINER_CONF_FILE={config}\nAPPTAINER_SUID_INSTALL={state["suid"]}\n'
        assert command[1:3] == ['inspect', '--json']
        assert Path(command[-1]).read_bytes() == b'complete SIF'
        state['inspections'].append(command[-1])
        return json.dumps({'data': {'attributes': {'labels': {}}}})
    def run(command, **kwargs):
        assert command[:2] == ['apptainer', 'pull']
        env = kwargs['env']
        assert env['APPTAINER_CONFIG_FILE'] == env['SINGULARITY_CONFIG_FILE']
        assert Path(env['APPTAINER_TMPDIR']).is_dir()
        state['environments'].append(dict(env))
        Path(env['APPTAINER_CACHEDIR']).mkdir(parents=True, exist_ok=True)
        (Path(env['APPTAINER_CACHEDIR']) / 'layer').write_bytes(b'OCI download')
        state['pulls'].append(Path(env['APPTAINER_CONFIG_FILE']).read_text())
        Path(command[-2]).write_bytes(b'partial' if state['fail'] else b'complete SIF')
        if state['fail']:
            raise subprocess.CalledProcessError(255, command)
        return subprocess.CompletedProcess(command, 0)
    monkeypatch.setattr(module.subprocess, 'check_output', output)
    monkeypatch.setattr(module.subprocess, 'run', run)
    monkeypatch.setattr(module.os, 'geteuid', lambda: 1001)
    return config, state


def test_bounded_pull_preserves_site_config_and_reuses_offline(tmp_path, runtime, monkeypatch):
    config, state = runtime
    before = config.read_bytes()
    image = tmp_path / 'cache/release.simg'
    supplied = {'TEST_PRESERVE': 'yes'}
    module.pull_image('apptainer', 'docker://release', image, cpus=2, env=supplied)
    assert image.read_bytes() == b'complete SIF'
    assert config.read_bytes() == before
    assert supplied == {'TEST_PRESERVE': 'yes'}
    assert 'mksquashfs procs = 2' in state['pulls'][0]
    assert 'mksquashfs mem = 3814M' in state['pulls'][0]
    assert len(state['inspections']) == 1
    monkeypatch.setattr(module, 'compression_config', lambda *a: pytest.fail('configuration on cache hit'))
    module.pull_image('apptainer', 'docker://release', image, env={})
    assert len(state['pulls']) == 1 and len(state['inspections']) == 2
    assert not list(image.parent.glob('.container-*'))


def test_failed_conversion_never_publishes_partial_cache(tmp_path, runtime):
    config, state = runtime
    image = tmp_path / 'cache/release.simg'
    state['fail'] = True
    with pytest.raises(subprocess.CalledProcessError):
        module.pull_image('apptainer', 'docker://release', image, env={})
    assert not image.exists() and not list(image.parent.glob('.container-*'))
    state['fail'] = False
    module.pull_image('apptainer', 'docker://release', image, env={})
    assert image.read_bytes() == b'complete SIF'


@pytest.mark.parametrize('suid', ['0', '1'])
def test_stricter_site_limits_are_preserved(tmp_path, runtime, suid):
    config, state = runtime
    config.write_text('mksquashfs procs = 1\nmksquashfs mem = 512M\n')
    state['suid'] = suid
    module.pull_image('apptainer', 'docker://release', tmp_path / 'image.sif', cpus=4, env={})
    assert 'mksquashfs procs = 1' in state['pulls'][0]
    assert 'mksquashfs mem = 512M' in state['pulls'][0]


def test_setuid_runtime_without_site_caps_stops_before_unbounded_pull(tmp_path, runtime):
    _, state = runtime
    state['suid'] = '1'
    with pytest.raises(ValueError, match='administrator mksquashfs'):
        module.pull_image('apptainer', 'docker://release', tmp_path / 'image.sif', env={})
    assert not state['pulls']


@pytest.mark.parametrize('fail', [False, True])
def test_image_expansion_and_default_layer_cache_stay_in_owned_task_scratch(tmp_path, runtime, fail):
    _, state = runtime
    state['fail'] = fail
    local = tmp_path / 'node'; local.mkdir()
    sentinel = local / 'other-job'; sentinel.write_text('keep')
    image = tmp_path / 'shared/image.sif'
    supplied = {'PHENORADAR_TASK_SCRATCH': str(local)}
    try:
        module.pull_image('apptainer', 'docker://release', image, env=supplied)
    except subprocess.CalledProcessError:
        assert fail
    env = state['environments'][0]
    assert Path(env['APPTAINER_TMPDIR']).is_relative_to(local)
    assert Path(env['APPTAINER_CACHEDIR']).is_relative_to(local)
    assert env['APPTAINER_CACHEDIR'] == env['SINGULARITY_CACHEDIR']
    assert list(local.iterdir()) == [sentinel]
    assert supplied == {'PHENORADAR_TASK_SCRATCH': str(local)}
    assert image.exists() is not fail


def test_explicit_layer_cache_is_preserved(tmp_path, runtime):
    _, state = runtime
    local, cache = tmp_path / 'node', tmp_path / 'custom-cache'
    local.mkdir(); cache.mkdir(); (cache / 'keep').write_text('caller owned')
    module.pull_image('apptainer', 'docker://release', tmp_path / 'image.sif',
                      env={'TMPDIR': str(local), 'APPTAINER_CACHEDIR': str(cache)})
    assert state['environments'][0]['APPTAINER_CACHEDIR'] == str(cache)
    assert (cache / 'keep').read_text() == 'caller owned' and (cache / 'layer').exists()
    assert not list(local.iterdir())
