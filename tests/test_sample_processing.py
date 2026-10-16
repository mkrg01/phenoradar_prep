"""Software upgrades keep completed samples and never mix native generating releases."""
from pathlib import Path
import shutil
import subprocess
import sys

import pytest
import yaml

import dataset
from build_dispatch import verify_and_stage
from build_products import load_complete
from common import read_tsv, write_tsv
from database_fixtures import database_from_stages
from dataset_assets import ProductIntegrityError
from sample_processing import FIELDS, generating_conditions, native_processing
from test_datasets import dataset_project, fake_genegalleon, imported, native_events


@pytest.fixture
def upgrade_project(dataset_project, monkeypatch):
    root = dataset_project
    (root / 'VERSION').write_text('0.4.9\n')
    repo = fake_genegalleon(root)
    (repo / 'VERSION').write_text('1.0.0\n')
    config = root / 'config/build.yaml'
    cfg = yaml.safe_load(config.read_text())
    cfg['genegalleon']['version'] = '1.0.0'
    config.write_text(yaml.safe_dump(cfg))
    old = dataset.prepare(root, 'old', config)
    dataset.submit(old, until='quant', dry_run=True)
    dataset.worker(old, 1)
    monkeypatch.setenv('FAKE_GG_FAIL_QUANT', '1')
    with pytest.raises(subprocess.CalledProcessError):
        dataset.worker(old, 2)
    monkeypatch.delenv('FAKE_GG_FAIL_QUANT')
    new_repo = root / 'new_gg'
    shutil.copytree(repo, new_repo)
    (new_repo / 'VERSION').write_text('1.0.1\n')
    (new_repo / 'genegalleon.sif').write_text('updated container')
    cfg['genegalleon'].update(version='1.0.1', revision='f' * 40, repository=str(new_repo))
    cfg['reuse_from'] = str(old)
    config.write_text(yaml.safe_dump(cfg))
    return root, config, old


@pytest.mark.parametrize('deferred', [False, True])
def test_upgrade_reuses_complete_sample_and_restarts_incomplete_sample(upgrade_project, deferred):
    root, config, old = upgrade_project
    report = dataset.plan(root, config, name='new')[-1]
    assert [row['assembly'] for row in report] == ['reuse', 'pending', 'pending']
    assert [row['quant'] for row in report] == ['reuse', 'pending', 'pending']
    build = dataset.prepare(root, 'new', config, defer=deferred)
    if deferred:
        verify_and_stage(build, dataset.load(build))
    manifest = dataset.load(build)
    old_manifest = dataset.load(old)
    first, second = manifest['items'][:2]
    origin = native_processing(dataset.item_products(old_manifest, old_manifest['items'][0]))
    assert native_processing(dataset.item_products(manifest, first)) == origin
    assert origin['genegalleon_version'].startswith('1.0.0+local.')
    assert not (build / 'work/genegalleon' / first['species']).exists()
    assert not (build / 'work/genegalleon' / second['species']).exists()
    commands = dataset.submit(build, until='quant', dry_run=True)
    assert '--array=2,3' in commands[0]
    dataset.worker(build, 1)
    dataset.worker(build, 2)
    assert len(native_events(build)) == 1
    event = native_events(build)[0]
    for flag in ('RUN_ASSEMBLY', 'RUN_BUSCO_LONGEST_CDS', 'RUN_AMALGKIT_QUANT'):
        assert event['env']['GG_TRANSCRIPTOME_' + flag] == '1'
    completed = native_processing(dataset.item_products(dataset.load(build), second))
    assert completed['genegalleon_version'].startswith('1.0.1+local.')
    assert len(native_events(old)) == 2


def test_same_software_still_resumes_registered_checkpoints(upgrade_project):
    root, config, old = upgrade_project
    cfg = yaml.safe_load(config.read_text())
    cfg['genegalleon'] = dataset.load(old)['config']['genegalleon']
    # Frozen resolution adds paths that were absent from the original config;
    # use its actual original local override and leave image selection implicit.
    cfg['genegalleon'].pop('image', None)
    config.write_text(yaml.safe_dump(cfg))
    report = dataset.plan(root, config, name='same')[-1]
    assert [report[1][s] for s in dataset.STAGES] == ['reuse', 'reuse', 'pending']


def test_lost_receipt_cannot_mix_an_old_completed_sample_with_new_software(upgrade_project):
    root, config, old = upgrade_project
    build = dataset.prepare(root, 'repair', config)
    manifest = dataset.load(build)
    item = manifest['items'][0]
    source = next((Path(manifest['config']['store']) / item['species']).glob('*/busco.json'))
    source.unlink()
    with pytest.raises(ProductIntegrityError, match='complete restart'):
        dataset.item_products(manifest, item)
    assert item['species'] in verify_and_stage(build, dataset.load(build))
    manifest = dataset.load(build)
    assert not any(dataset.item_products(manifest, manifest['items'][0]).values())
    dataset.submit(build, until='quant', dry_run=True)
    dataset.worker(build, 1)
    assert native_events(build)[0]['env']['GG_TRANSCRIPTOME_RUN_ASSEMBLY'] == '1'
    assert native_processing(dataset.item_products(dataset.load(build), item))[
        'genegalleon_version'].startswith('1.0.1+local.')


def test_completed_database_preserves_origin_through_multiple_reuses(upgrade_project):
    root, config, old = upgrade_project
    manifest = dataset.load(old)
    source = database_from_stages(root, manifest['config']['store'], manifest['items'][:1])
    original = read_tsv(source / 'sample_processing.tsv')[0]
    assert original['genegalleon_version'].startswith('1.0.0+local.')
    cfg = yaml.safe_load(config.read_text())
    cfg['reuse_from'] = str(source)
    config.write_text(yaml.safe_dump(cfg))
    item = manifest['items'][0]
    metadata = root / 'input/one.tsv'
    write_tsv(metadata, manifest['fields'], [item['row']])
    first = dataset.prepare(root, 'reuse_once', config, metadata)
    first_manifest = dataset.load(first)
    processing = native_processing(dataset.item_products(first_manifest, first_manifest['items'][0]))
    assert processing == load_complete(source)['products'][item['species']]['processing']
    second_source = database_from_stages(root, first_manifest['config']['store'],
                                         first_manifest['items'], 'second_source')
    cfg['reuse_from'] = str(second_source)
    config.write_text(yaml.safe_dump(cfg))
    second = dataset.prepare(root, 'reuse_twice', config, metadata, defer=True)
    verify_and_stage(second, dataset.load(second))
    second_manifest = dataset.load(second)
    assert native_processing(dataset.item_products(second_manifest, second_manifest['items'][0])) == processing
    assert read_tsv(second_source / 'sample_processing.tsv')[0]['genegalleon_version'] == original['genegalleon_version']


def test_explicit_sample_rebuild_skips_reuse_without_touching_source(upgrade_project):
    root, config, old = upgrade_project
    dataset.worker(old, 2)
    dataset.worker(old, 3)
    before = (old / 'build.json').read_bytes()
    cfg = yaml.safe_load(config.read_text())
    cfg['rebuild_samples'] = ['A1']
    config.write_text(yaml.safe_dump(cfg))
    report = dataset.plan(root, config, name='forced')[-1]
    assert [r['assembly'] for r in report] == ['pending', 'reuse', 'reuse']
    build = dataset.prepare(root, 'forced', config, defer=True)
    verify_and_stage(build, dataset.load(build))
    assert [r['assembly'] for r in dataset.status(build)] == ['pending', 'reuse', 'reuse']
    dataset.submit(build, until='quant', dry_run=True)
    dataset.worker(build, 1)
    assert native_events(build)[0]['env']['GG_TRANSCRIPTOME_RUN_ASSEMBLY'] == '1'
    assert (old / 'build.json').read_bytes() == before
    cfg.update(reuse_from=None, rebuild_samples=[])
    config.write_text(yaml.safe_dump(cfg))
    assert all(r['assembly'] == 'pending' for r in dataset.plan(root, config, name='fresh')[-1])


@pytest.mark.parametrize('selection', ['A1', ['missing'], [1]])
def test_invalid_rebuild_selection_is_rejected(dataset_project, selection):
    config = dataset_project / 'config/build.yaml'
    cfg = yaml.safe_load(config.read_text())
    cfg['rebuild_samples'] = selection
    config.write_text(yaml.safe_dump(cfg))
    with pytest.raises(ValueError, match='rebuild_samples'):
        dataset.plan(dataset_project, config)


def test_reporting_and_release_labels_do_not_block_resume(upgrade_project):
    root, _, old = upgrade_project
    reporting = root / 'workflow/scripts/build_reports.py'
    reporting.write_text(reporting.read_text() + '\n# Presentation-only update.\n')
    (root / 'VERSION').write_text('0.4.10\n')
    dataset.load(old, check_code=True)
    assert dataset.submit(old, until='quant', dry_run=True)
    scientific = root / 'workflow/scripts/aggregate_tpm.py'
    scientific.write_text(scientific.read_text() + '\n# Processing code changed.\n')
    with pytest.raises(ProductIntegrityError, match='registered file changed'):
        dataset.load(old, check_code=True)


def test_mixed_generating_conditions_cannot_be_reported_as_one_release():
    products = {key: {'provenance': {'processing': {'genegalleon_version': version}}}
                for key, version in zip(('reference', 'busco', 'quant'), ('1.0.0', '1.0.0', '1.0.1'))}
    with pytest.raises(ValueError, match='different GeneGalleon'):
        native_processing(products)


def test_local_overrides_do_not_claim_requested_release_pins(dataset_project):
    cfg = yaml.safe_load((dataset_project / 'config/build.yaml').read_text())
    cfg['genegalleon']['repository'] = str(dataset_project / 'local')
    lock = {'source': {'kind': 'local_source', 'version': '1.0.0',
                       'files': [{'path': str(dataset_project / 'local/VERSION'), 'sha256': 'local-source'}]},
            'container': {'kind': 'local_image', 'files': [{'sha256': 'local-image'}]}}
    processing = generating_conditions({'config': cfg, 'software_lock': lock})
    assert processing['genegalleon_version'].startswith('1.0.0+local.')
    assert processing['revision'] is None
    assert processing['image_uri'] is None
    assert processing['container_sha256'] == ['local-image']


@pytest.mark.parametrize('reverse', [False, True])
def test_complete_sample_wins_over_partial_work_from_another_release(reverse):
    from build_stage_reuse import merge
    complete = {'kind': 'stages', 'processing': {'genegalleon_version': '1.0.0'},
                'products': {'reference': {'cds': {'sha256': 'old-cds'}},
                             'busco': {'full': {'sha256': 'busco'}},
                             'quant': {'abundance': {'sha256': 'abundance'}}}}
    partial = {'kind': 'stages', 'processing': {'genegalleon_version': '1.0.1'},
               'products': {'reference': {'cds': {'sha256': 'new-cds'}}, 'busco': None, 'quant': None}}
    selected = {}
    for entry in ([partial, complete] if reverse else [complete, partial]):
        merge(selected, 'sample', entry)
    assert selected['sample'] is complete


def test_legacy_export_is_standalone_and_never_invents_a_version(dataset_project, monkeypatch, capsys):
    root = dataset_project
    source = imported(root)
    before = (source / 'manifest.json').read_bytes()
    output = root / 'supplementary'
    monkeypatch.setattr(sys, 'argv', ['dataset.py', 'export-processing', '--build', str(source),
                                      '--outdir', str(output)])
    assert dataset.main() == 0
    rows = read_tsv(output / 'sample_processing.tsv')
    assert len(rows) == 3
    assert list(rows[0]) == FIELDS
    assert {row['genegalleon_version'] for row in rows} == {'unknown'}
    assert {row['assembly_method'] for row in rows} == {'unknown'}
    assert rows[0]['busco_complete_percent'] == '80.00'
    assert rows[0]['orthodb_clade'] == 'Embryophyta'
    assert (source / 'manifest.json').read_bytes() == before
    assert not (root / 'results').exists()
    (source / 'sample_processing.tsv').write_text('changed summary')
    with pytest.raises(ProductIntegrityError, match='registered file changed'):
        load_complete(source)


def test_scientific_changes_still_prevent_cross_version_reuse(upgrade_project):
    root, config, old = upgrade_project
    manifest = dataset.load(old)
    source = database_from_stages(root, manifest['config']['store'], manifest['items'][:1])
    cfg = yaml.safe_load(config.read_text())
    cfg['reuse_from'] = str(source)
    cfg['genegalleon']['settings']['assembly_method'] = 'Trinity'
    config.write_text(yaml.safe_dump(cfg))
    assert dataset.plan(root, config, name='different')[-1][0]['assembly'] == 'conflict'
