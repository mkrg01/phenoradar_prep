"""Shared references publish atomically, stay pinned, and reach native jobs."""
from concurrent.futures import ThreadPoolExecutor
import csv
import json
from pathlib import Path
import subprocess
import threading

import pytest
import yaml

import busco_reference as busco
from build_dispatch import verify_and_stage
from dataset import load, prepare, submit, worker
from dataset_assets import ProductIntegrityError, UncachedVerification
from test_datasets import dataset_project, fake_genegalleon, native_events, new_dataset


@pytest.fixture
def container(monkeypatch):
    calls = {'downloads': [], 'probes': [], 'native': [], 'fail_download': False, 'fail_mount': False}
    run = subprocess.run
    monkeypatch.setattr(busco, 'runtime', lambda: '/fake/apptainer')

    def execute(command, **kwargs):
        if command[0] != '/fake/apptainer':
            if str(command[-1]).endswith('gg_transcriptome_generation_entrypoint.sh'):
                calls['native'].append(kwargs['env'])
            return run(command, **kwargs)
        if '--mount' in command:
            spec = next(csv.reader([command[command.index('--mount') + 1]]))
            mount = dict(field.split('=', 1) for field in spec)
            staging = Path(mount['src'])
            lineage = command[-1]
            assert not (staging.parent / lineage).exists()  # No partial public snapshot.
            calls['downloads'].append(staging)
            out = staging / 'dataset/lineages' / lineage
            (out / 'hmms').mkdir(parents=True)
            (out / 'dataset.cfg').write_text('[busco_dataset]\ncreation_date=2026-01-01\n')
            (out / 'hmms/marker.hmm').write_text('synthetic HMM\n')
            if calls['fail_download']:
                raise subprocess.CalledProcessError(1, command)
        else:
            calls['probes'].append(kwargs['env'])
            assert command[-2].startswith(busco.DESTINATION + '/lineages/')
            if calls['fail_mount']:
                raise subprocess.CalledProcessError(1, command)
        return subprocess.CompletedProcess(command, 0)

    monkeypatch.setattr(busco.subprocess, 'run', execute)
    return calls


def config(tmp_path):
    return {'lineage': 'embryophyta_odb12', 'cache_dir': str(tmp_path / 'busco')}


def test_parallel_preparation_downloads_once_and_reuses_offline(tmp_path, container):
    cfg = config(tmp_path)
    start = threading.Barrier(2)

    def acquire():
        start.wait(timeout=10)
        return busco.prepare(cfg, {'image': 'pinned.sif'})

    with ThreadPoolExecutor(max_workers=2) as pool:
        refs = list(pool.map(lambda _: acquire(), range(2)))
    assert refs[0] == refs[1]
    assert len(container['downloads']) == 1
    assert refs[0]['dataset_version'] == '2026-01-01'
    assert busco.prepare(cfg, {'image': 'unavailable.sif'}, expected=refs[0]) == refs[0]
    assert len(container['downloads']) == 1


def test_failed_download_is_never_published_and_retry_succeeds(tmp_path, container):
    cfg = config(tmp_path)
    container['fail_download'] = True
    with pytest.raises(subprocess.CalledProcessError):
        busco.prepare(cfg, {'image': 'pinned.sif'})
    assert not (Path(cfg['cache_dir']) / cfg['lineage']).exists()
    assert not any(Path(cfg['cache_dir']).glob('.*'))
    container['fail_download'] = False
    busco.prepare(cfg, {'image': 'pinned.sif'})
    assert len(container['downloads']) == 2


@pytest.mark.parametrize('change', ['content', 'missing', 'extra', 'receipt'])
def test_reference_damage_fails_instead_of_redownloading(tmp_path, container, change):
    cfg = config(tmp_path)
    ref = busco.prepare(cfg, {'image': 'pinned.sif'})
    marker = Path(ref['dataset']) / 'lineages/embryophyta_odb12/hmms/marker.hmm'
    if change == 'content':
        marker.write_text('corrupted HMM\n')
    elif change == 'missing':
        marker.unlink()
    elif change == 'extra':
        marker.with_name('unexpected.hmm').write_text('extra')
    else:
        Path(ref['receipt']['path']).write_text('{}')
    with pytest.raises(ProductIntegrityError):
        busco.prepare(cfg, {'image': 'pinned.sif'}, expected=ref)
    assert len(container['downloads']) == 1


def test_cached_only_does_not_download_or_create_directories(tmp_path, container):
    cfg = config(tmp_path)
    with pytest.raises(UncachedVerification):
        busco.prepare(cfg, {'image': 'pinned.sif'}, cached_only=True)
    assert not Path(cfg['cache_dir']).exists()
    assert container['downloads'] == []


def test_mount_preserves_other_mounts_and_quotes_paths(container):
    env = {'APPTAINER_MOUNT': 'type=bind,src=/a,dst=/data',
           'SINGULARITY_MOUNT': 'type=bind,src=/b,dst=/other,ro'}
    busco.mount_environment(env, {'dataset': '/shared/path, with spaces'})
    mounts = list(csv.reader(env['APPTAINER_MOUNT'].splitlines()))
    assert len(mounts) == 3
    assert mounts[-1] == ['type=bind', 'src=/shared/path, with spaces', 'dst=' + busco.DESTINATION, 'ro']
    assert 'SINGULARITY_MOUNT' not in env
    with pytest.raises(ValueError, match='overlaps'):
        busco.mount_environment({'APPTAINER_MOUNT': 'type=bind,src=/a,dst=/usr/local'}, {'dataset': '/b'})


@pytest.mark.parametrize('deferred', [False, True])
def test_build_prepares_once_and_mounts_same_reference_for_each_sample(dataset_project, container, deferred):
    root = dataset_project
    cfg_path = root / 'config/build.yaml'
    cfg = yaml.safe_load(cfg_path.read_text())
    cfg['busco']['cache_dir'] = 'resources/busco'
    cfg_path.write_text(yaml.safe_dump(cfg))
    if deferred:
        fake_genegalleon(root)
        build = prepare(root, 'shared', cfg_path, defer=True)
        assert not container['downloads']
        verify_and_stage(build, load(build))
    else:
        build = new_dataset(root, ('New plant', 'Other plant'))
    manifest = load(build)
    ref = manifest['busco_reference']
    assert len(container['downloads']) == 1
    submit(build, until='quant', dry_run=True)
    for index in (1, 2):
        worker(build, index)
    assert len(container['downloads']) == 1
    assert len(container['probes']) == len(container['native']) == 2
    for env in container['native']:
        assert 'src=' + ref['dataset'] in env['APPTAINER_MOUNT']
        assert env['APPTAINER_MOUNT'].endswith(',ro')
    from dataset import item_products
    products = item_products(load(build), manifest['items'][0])
    assert products['busco']['provenance']['busco_reference']['sha256'] == ref['sha256']
    assert not list((build / 'work/genegalleon').glob('*/downloads/busco_downloads'))


def test_mount_failure_stops_before_native_job(dataset_project, container):
    root = dataset_project
    cfg_path = root / 'config/build.yaml'
    cfg = yaml.safe_load(cfg_path.read_text())
    cfg['busco']['cache_dir'] = 'resources/busco'
    cfg_path.write_text(yaml.safe_dump(cfg))
    build = new_dataset(root)
    submit(build, until='quant', dry_run=True)
    container['fail_mount'] = True
    with pytest.raises(subprocess.CalledProcessError):
        worker(build, 1)
    assert not native_events(build)
    assert not container['native']
    receipt = json.loads((build / 'jobs/status/New_plant_SRR1.sample.json').read_text())
    assert receipt['state'] == 'failed'


@pytest.mark.parametrize('value', ['../escape', '/local/lineage', '', None])
def test_lineage_requires_a_dataset_name(tmp_path, value):
    with pytest.raises(ValueError, match='dataset name'):
        busco.settings(tmp_path, {'lineage': value})


@pytest.mark.parametrize('version,help_text', [('singularity-ce version 4.0', '--mount'),
                                              ('apptainer version 1.0', '--bind')])
def test_unsupported_runtime_fails_clearly(monkeypatch, version, help_text):
    busco.runtime.cache_clear()
    monkeypatch.setattr(busco.shutil, 'which', lambda name: '/fake/' + name)
    monkeypatch.setattr(busco.subprocess, 'check_output',
                        lambda command, **kwargs: version if command[-1] == '--version' else help_text)
    with pytest.raises(ValueError, match='APPTAINER_MOUNT support'):
        busco.runtime()
    busco.runtime.cache_clear()


def test_nested_mount_cannot_shadow_lineage_files(container):
    with pytest.raises(ValueError, match='overlaps'):
        busco.mount_environment({'APPTAINER_MOUNT': 'type=bind,src=/a,target=' + busco.DESTINATION + '/lineages'},
                                {'dataset': '/shared'})
