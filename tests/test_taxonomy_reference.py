"""Taxonomy snapshots are shared, immutable, and checked before native jobs."""
from concurrent.futures import ThreadPoolExecutor
import csv
import json
from pathlib import Path
import subprocess
import threading

import pytest
import yaml

import busco_reference as busco
import taxonomy_reference as taxonomy
from build_dispatch import verify_and_stage
from dataset import item_products, load, prepare, submit, worker
from dataset_assets import ProductIntegrityError, UncachedVerification
from test_datasets import dataset_project, fake_genegalleon, native_events, new_dataset


@pytest.fixture
def container(monkeypatch):
    calls = {'taxonomy': 0, 'busco': 0, 'probes': [], 'native': [], 'fail_download': False, 'fail_mount': False}
    run = subprocess.run
    monkeypatch.setattr(busco, 'runtime', lambda: '/fake/apptainer')

    def execute(command, **kwargs):
        if command[0] != '/fake/apptainer':
            if str(command[-1]).endswith('gg_transcriptome_generation_entrypoint.sh'):
                calls['native'].append(kwargs['env'])
            return run(command, **kwargs)
        mounts = [dict(field.split('=', 1) for field in next(csv.reader([command[i+1]])) if '=' in field)
                  for i, value in enumerate(command) if value == '--mount']
        if mounts and mounts[0]['dst'] == '/taxonomy':
            calls['taxonomy'] += 1
            staging = Path(mounts[0]['src'])
            assert not (staging.parent / 'ete').exists()
            directory = staging / 'work/downloads/ete_taxonomy'
            for name in ('ete', 'ete4'):
                (directory / name).mkdir(parents=True)
            for name in taxonomy.REQUIRED:
                (directory / name).write_text('synthetic ' + name)
            if calls['fail_download']:
                raise subprocess.CalledProcessError(1, command)
        elif mounts and mounts[0]['dst'] == '/busco':
            calls['busco'] += 1
            directory = Path(mounts[0]['src']) / 'dataset/lineages' / command[-1]
            (directory / 'hmms').mkdir(parents=True)
            (directory / 'dataset.cfg').write_text('creation_date=2026-01-01\n')
            (directory / 'hmms/marker.hmm').write_text('synthetic HMM')
        else:
            calls['probes'].append(command)
            if calls['fail_mount']:
                raise subprocess.CalledProcessError(1, command)
        return subprocess.CompletedProcess(command, 0)

    monkeypatch.setattr(subprocess, 'run', execute)
    return calls


def config(tmp_path):
    return {'cache_dir': str(tmp_path / 'taxonomy')}


def test_parallel_preparation_publishes_once_and_reuses_offline(tmp_path, container):
    cfg = config(tmp_path)
    software = {'repository': '/pinned/source', 'image': 'pinned.sif'}
    start = threading.Barrier(2)

    def acquire(_):
        start.wait(timeout=10)
        return taxonomy.prepare(cfg, software)

    with ThreadPoolExecutor(max_workers=2) as pool:
        refs = list(pool.map(acquire, range(2)))
    assert refs[0] == refs[1] and container['taxonomy'] == 1
    assert taxonomy.prepare(cfg, {}, expected=refs[0], cached_only=True) == refs[0]
    assert container['taxonomy'] == 1


def test_failed_preparation_never_publishes_partial_data(tmp_path, container):
    cfg = config(tmp_path)
    container['fail_download'] = True
    with pytest.raises(subprocess.CalledProcessError):
        taxonomy.prepare(cfg, {'repository': '/pinned/source', 'image': 'pinned.sif'})
    assert not (Path(cfg['cache_dir']) / 'ete').exists()
    assert not list(Path(cfg['cache_dir']).glob('.ete-*'))
    container['fail_download'] = False
    taxonomy.prepare(cfg, {'repository': '/pinned/source', 'image': 'pinned.sif'})
    assert container['taxonomy'] == 2


@pytest.mark.parametrize('change', ['content', 'missing', 'extra', 'symlink', 'receipt', 'directory'])
def test_damage_stops_without_refreshing_snapshot(tmp_path, container, change):
    cfg = config(tmp_path)
    ref = taxonomy.prepare(cfg, {'repository': '/pinned/source', 'image': 'pinned.sif'})
    target = Path(ref['dataset']) / 'taxa.sqlite'
    if change == 'content':
        target.write_text('changed')
    elif change == 'missing':
        target.unlink()
    elif change == 'extra':
        target.with_name('unexpected').write_text('extra')
    elif change == 'symlink':
        original = tmp_path / 'external.sqlite'
        target.rename(original)
        target.symlink_to(original)
    elif change == 'directory':
        (target.parent / 'ete4').rmdir()
    else:
        Path(ref['receipt']['path']).write_text('{}')
    with pytest.raises(ProductIntegrityError):
        taxonomy.prepare(cfg, {}, expected=ref)
    assert container['taxonomy'] == 1


def test_cached_only_missing_snapshot_creates_nothing(tmp_path, container):
    cfg = config(tmp_path)
    with pytest.raises(UncachedVerification):
        taxonomy.prepare(cfg, {}, cached_only=True)
    assert not Path(cfg['cache_dir']).exists()


def test_disabled_and_legacy_configs_need_no_reference(container):
    assert taxonomy.prepare({'cache_dir': None}, {}) is None
    assert taxonomy.prepare({}, {}) is None
    assert container['taxonomy'] == 0


@pytest.mark.parametrize('deferred', [False, True])
def test_build_shares_busco_and_taxonomy_and_records_provenance(dataset_project, container, deferred):
    root = dataset_project
    config_path = root / 'config/build.yaml'
    cfg = yaml.safe_load(config_path.read_text())
    cfg['taxonomy']['cache_dir'] = 'resources/taxonomy'
    cfg['busco']['cache_dir'] = 'resources/busco'
    config_path.write_text(yaml.safe_dump(cfg))
    if deferred:
        fake_genegalleon(root)
        build = prepare(root, 'shared', config_path, defer=True)
        assert container['taxonomy'] == 0
        verify_and_stage(build, load(build))
    else:
        build = new_dataset(root, ('New plant', 'Other plant'))
    manifest = load(build)
    submit(build, until='quant', dry_run=True)
    for index in range(1, len(manifest['items']) + 1):
        worker(build, index)
    assert container['taxonomy'] == container['busco'] == 1
    assert len(container['native']) == len(manifest['items'])
    assert len(container['probes']) == 2 * len(manifest['items'])
    for env in container['native']:
        mounts = list(csv.reader(env['APPTAINER_MOUNT'].splitlines()))
        assert {taxonomy.DESTINATION, busco.DESTINATION} == {
            field[4:] for row in mounts for field in row if field.startswith('dst=')}
        assert all(row[-1] == 'ro' for row in mounts)
    for item in manifest['items']:
        products = item_products(load(build), item)
        for key in ('reference', 'busco', 'quant'):
            assert products[key]['provenance']['taxonomy_reference'] == manifest['taxonomy_reference']
    assert not list((build / 'work/genegalleon').glob('*/downloads/ete_taxonomy'))


def test_bad_mount_records_failure_without_starting_native_job(dataset_project, container):
    path = dataset_project / 'config/build.yaml'
    cfg = yaml.safe_load(path.read_text())
    cfg['taxonomy']['cache_dir'] = 'resources/taxonomy'
    path.write_text(yaml.safe_dump(cfg))
    build = new_dataset(dataset_project)
    submit(build, until='quant', dry_run=True)
    container['fail_mount'] = True
    with pytest.raises(subprocess.CalledProcessError):
        worker(build, 1)
    assert not native_events(build)
    assert json.loads((build / 'jobs/status/New_plant_SRR1.sample.json').read_text())['state'] == 'failed'


@pytest.mark.parametrize('value', ['', True, '\nunsafe'])
def test_invalid_cache_paths_fail(tmp_path, value):
    with pytest.raises(ValueError, match='taxonomy.cache_dir'):
        taxonomy.settings(tmp_path, {'cache_dir': value})


def test_mount_cannot_be_shadowed_and_preserves_unrelated_mounts(container):
    env = {'SINGULARITY_MOUNT': 'type=bind,src=/other,dst=/data'}
    taxonomy.mount_environment(env, {'dataset': '/shared/path, with spaces'})
    assert len(list(csv.reader(env['APPTAINER_MOUNT'].splitlines()))) == 2
    with pytest.raises(ValueError, match='overlaps'):
        taxonomy.mount_environment({'APPTAINER_MOUNT': 'type=bind,src=/bad,dst=/workspace'}, {'dataset': '/shared'})
