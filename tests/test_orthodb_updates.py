"""Refresh one shared OrthoDB reference without repeating native sample processing."""
import json
import os
from pathlib import Path
import shutil

import pytest
import yaml

import dataset
from build_dispatch import verify_and_stage
from build_products import complete, load_complete
from common import read_tsv, sha256, write_json, write_tsv
from configuration import validate_keys
from incremental_odb import plan as plan_mapping
from mapping_tables import load_tables, reference_for
from odb_settings import pin_reference, reference_directory
from prepare_odb_reference import prepare as prepare_reference
from run_odb_chunk import run as map_chunk
from test_datasets import dataset_project, imported
from test_portable_build import execute
from test_incremental_odb import snapshot
from database_fixtures import database_from_stages
from dataset_assets import identities


def configure(root, **odb):
    config = root / 'config/build.yaml'
    cfg = yaml.safe_load(config.read_text())
    cfg['odb'].update(odb)
    config.write_text(yaml.safe_dump(cfg))
    return config


@pytest.mark.parametrize('settings', [
    {'version': 'current'}, {'version': 'latest'}, {'version': None},
    {'version': 'v12/other'}, {'snapshot': None}, {'snapshot': '2026-10-16'},
    {'force_remap': False}, {'force_remap': True},
])
def test_reject_moving_releases_and_invalid_settings(settings):
    with pytest.raises(ValueError, match='odb\\.'):
        validate_keys({'odb': settings})


@pytest.mark.parametrize('change', [
    {'version': 'v11'}, {'ncbi_tax_id': 33090},
])
@pytest.mark.parametrize('deferred', [False, True])
def test_refresh_imports_native_and_translation_without_old_mapping(dataset_project, change, deferred, monkeypatch):
    root = dataset_project
    source = imported(root)
    old = load_complete(source)
    manifest_bytes = (source / 'manifest.json').read_bytes()
    # Remapping must not require old mapping or expression payloads, even when
    # only the portable database, without its original working tree, is available.
    Path(old['mapping']['path']).unlink()
    for product in old['products'].values():
        for key in ('expression', 'expression_qc'):
            Path(product[key]['path']).unlink()
    config = configure(root, **change)
    def forbidden(*args, **kwargs):
        pytest.fail('OrthoDB refresh must not prepare GeneGalleon software')
    monkeypatch.setattr(dataset, 'resolve_software', forbidden)
    report = dataset.plan(root, config)[-1]
    assert all(row[stage] == 'reuse' for row in report for stage in dataset.STAGES)
    assert all(row['mapping'] == 'pending_inputs' for row in report)
    build = dataset.prepare(root, 'refresh', config, defer=deferred)
    if deferred:
        verify_and_stage(build, dataset.load(build))
    manifest = dataset.load(build)
    for item in manifest['items']:
        products = dataset.item_products(manifest, item)
        assert all(products.values())
        assert products['reference']['cds']['sha256'] == old['products'][item['species']]['cds']['sha256']
    assert len(list((build / 'work/cache/proteins').glob('*/table_1/receipt.json'))) == len(report)
    assert not list((build / 'work/cache/odb').rglob('snapshot.json'))
    assert not list((build / 'work/cache/expression').rglob('receipt.json'))
    assert (source / 'manifest.json').read_bytes() == manifest_bytes


def test_release_changes_remap_and_same_release_reuses_results(dataset_project, fake_odb, command_environment, monkeypatch):
    root = dataset_project
    if not shutil.which('snakemake') or not shutil.which('seqkit'):
        pytest.skip('Snakemake and seqkit required')
    source = imported(root)
    original = load_complete(source)
    before = (source / 'manifest.json').read_bytes()
    taxonomy = root / 'resources/taxonomy/taxa.sqlite'
    taxonomy.parent.mkdir(parents=True)
    shutil.copy2(root / 'input/taxa.sqlite', taxonomy)
    env = command_environment({'python': os.sys.executable, 'seqkit': shutil.which('seqkit'), 'ODB-mapper': fake_odb})
    events = root / 'mapping_events.txt'
    env['FAKE_ODB_LOG'] = str(events)
    def forbidden(*args, **kwargs):
        pytest.fail('complete native products must not run GeneGalleon')
    monkeypatch.setattr(dataset, 'resolve_software', forbidden)
    sample_count = len(original['products'])
    for index, (version, batches) in enumerate([
        ('v11', 1), ('v12', 2), ('v12', 2),
    ]):
        config = configure(root, version=version, chunk_size=1)
        build = dataset.prepare(root, f'refresh_{index}', config)
        dataset.materialize(build)
        execute(root, build, 'database', env, mapping=True)
        # Replanning an interrupted build retains its completed chunks.
        if index == 0:
            shutil.rmtree(build / 'work/database/orthogroups/mapping/incremental_plan')
            (build / 'work/database/orthogroups/mapping/snapshot.json').unlink()
        execute(root, build, 'database', env, mapping=True)
        complete(build)
        result = load_complete(build)
        if index == 0:
            import analysis
            probe = root / 'analysis-probe.yaml'
            probe.write_text(yaml.safe_dump({'inputs': {'species_trait': None},
                'phylogeny': {'trees': [], 'contrast_pairs': {'enabled': False}}}))
            pipeline = analysis.settings(root, probe, build)[1]
            assert pipeline['odb']['version'] == version
            assert 'snapshot' not in pipeline['odb']
        assert len(events.read_text().splitlines()) == batches * sample_count
        assert result['odb']['version'] == version
        assert 'snapshot' not in result['odb']
        mapping = json.loads(Path(result['mapping']['path']).read_text())
        assert mapping['reference_sha256s'] == [result['odb']['reference_sha256']]
        pin = json.loads((build / 'database/provenance/odb_reference_pin.json').read_text())
        if index < 2:
            assert pin['reference_sha256'] == result['odb']['reference_sha256']
        else:
            assert pin['reference_sha256'] is None  # No new mapping jobs in this build.
        for name, product in result['products'].items():
            for key in ('cds', 'busco', 'abundance', 'protein'):
                assert product[key]['sha256'] == original['products'][name][key]['sha256']
            translated = json.loads((build / f'work/database/proteins/{product["odb_species"]}_protein.json').read_text())
            assert translated['reused']
        summary = read_tsv(build / 'database/sample_processing.tsv')
        assert {r['orthodb_release'] for r in summary} == {version}
        assert all('orthodb_snapshot' not in row for row in summary)
        cfg = yaml.safe_load(config.read_text())
        cfg['reuse_from'] = str(build / 'database')
        config.write_text(yaml.safe_dump(cfg))
    assert (source / 'manifest.json').read_bytes() == before


@pytest.mark.parametrize('reverse', [False, True])
def test_mixed_release_sources_reuse_matching_samples_and_prefer_matching_duplicates(dataset_project, reverse):
    root = dataset_project
    imported(root)
    items = identities(root / 'input/metadata.tsv')[1]
    store = root / 'resources/dataset_assets'
    old = database_from_stages(root, store, items, name='old_release', version='v11')
    current = database_from_stages(root, store, items[:1], name='current_release')
    config = configure(root, version='v12')
    cfg = yaml.safe_load(config.read_text())
    cfg['reuse_from'] = [str(path) for path in ([current, old] if reverse else [old, current])]
    config.write_text(yaml.safe_dump(cfg))
    report = dataset.plan(root, config)[-1]
    assert all(row[stage] == 'reuse' for row in report for stage in dataset.STAGES)
    assert report[0]['mapping'] == 'reuse'
    assert all(row['mapping'] == 'pending_inputs' for row in report[1:])
    build = dataset.prepare(root, 'mixed_release', config, defer=True)
    verify_and_stage(build, dataset.load(build))
    cached = load_tables(build / 'work/cache/odb/v12_3193/reused/snapshot.json')
    assert set(cached['tables']) == {items[0]['species']}
    assert len(list((build / 'work/cache/expression').glob('*/receipt.json'))) == 1


def test_mixed_release_build_maps_only_outdated_samples(dataset_project, fake_odb, command_environment):
    root = dataset_project
    if not shutil.which('snakemake') or not shutil.which('seqkit'):
        pytest.skip('Snakemake and seqkit required')
    imported(root)
    items = identities(root / 'input/metadata.tsv')[1]
    store = root / 'resources/dataset_assets'
    current = database_from_stages(root, store, items[:1], name='current_release')
    old = database_from_stages(root, store, items[1:], name='old_release', version='v11')
    config = configure(root, version='v12', chunk_size=1)
    cfg = yaml.safe_load(config.read_text())
    cfg['reuse_from'] = [str(current), str(old)]
    config.write_text(yaml.safe_dump(cfg))
    taxonomy = root / 'resources/taxonomy/taxa.sqlite'
    taxonomy.parent.mkdir(parents=True)
    shutil.copy2(root / 'input/taxa.sqlite', taxonomy)
    env = command_environment({'python': os.sys.executable, 'seqkit': shutil.which('seqkit'), 'ODB-mapper': fake_odb})
    events = root / 'mapping_events.txt'
    env['FAKE_ODB_LOG'] = str(events)
    build = dataset.prepare(root, 'mixed_release', config)
    dataset.materialize(build)
    execute(root, build, 'database', env, mapping=True)
    complete(build)
    result = load_complete(build)
    mapping = load_tables(result['mapping']['path'])
    plan = json.loads((build / 'work/database/orthogroups/mapping/incremental_plan/plan.json').read_text())
    assert plan['reused_species'] == [items[0]['species']]
    assert plan['mapped_species'] == sorted(item['species'] for item in items[1:])
    assert len(events.read_text().splitlines()) == len(items) - 1
    assert mapping['version'] == 'v12'
    assert {row['orthodb_release'] for row in read_tsv(build / 'database/sample_processing.tsv')} == {'v12'}
    # Historical unknown checksums must not be inferred from other samples.
    assert reference_for(mapping, items[0]['species']) is None
    assert all(reference_for(mapping, item['species']) for item in items[1:])
    assert 'reference_sha256' not in result['odb']
    assert len(mapping['reference_sha256s']) == 1
    # A subsequent build of the same release imports every completed mapping.
    cfg['reuse_from'] = str(build / 'database')
    config.write_text(yaml.safe_dump(cfg))
    assert all(row['mapping'] == 'reuse' for row in dataset.plan(root, config)[-1])


def test_mapper_fallback_is_rejected_before_download(tmp_path, fake_odb, monkeypatch):
    monkeypatch.setenv('FAKE_ODB_API_VERSION', 'v12.0')
    with pytest.raises(ValueError, match='requested v11'):
        prepare_reference(tmp_path / 'ref', command=str(fake_odb), version='v11')
    assert not (tmp_path / 'ref/reference.json').exists()
    assert not list((tmp_path / 'ref').rglob('tiny.db'))


def test_shared_pin_blocks_reference_replacement_and_changed_release(tmp_path, fake_odb, frozen_reference, monkeypatch):
    marker = frozen_reference / 'reference.json'
    pin = tmp_path / 'pin.json'
    expected = sha256(marker)
    pin_reference(pin, 'v12', 3193, expected)
    protein = tmp_path / 'Alpha_plant_protein.fa'
    protein.write_text('>Alpha_plant_g1\nMK*\n')
    manifest = tmp_path / 'chunk_000.fs'
    manifest.write_text(str(protein) + '\n')
    value = json.loads(marker.read_text())
    value['changed'] = True
    write_json(marker, value)
    with pytest.raises(ValueError, match='pinned build snapshot'):
        map_chunk(manifest, marker, tmp_path / 'out', tmp_path / 'work', 'chunk_000',
                  command=str(fake_odb), jobs=1, batch_size=1, reference_pin=pin)
    assert not (tmp_path / 'out').exists()
    with pytest.raises(ValueError, match='version/node changed'):
        pin_reference(pin, 'v11', 3193)
    assert json.loads(pin.read_text())['reference_sha256'] == expected


def test_same_release_historical_mapping_is_reused_when_adding_samples(tmp_path):
    proteins = tmp_path / 'proteins'
    proteins.mkdir()
    for species in ('Alpha_plant', 'Beta_plant'):
        (proteins / f'{species}_protein.fa').write_text(f'>{species}_g1\nMK\n')
    snapshot(tmp_path / 'cache/v12_3193/old', proteins, ['Alpha_plant'])
    samples = tmp_path / 'samples.tsv'
    write_tsv(samples, ['species', 'odb_species'],
              [{'species': s, 'odb_species': s} for s in ('Alpha_plant', 'Beta_plant')])
    marker = tmp_path / 'reference.json'
    write_json(marker, {'version': 'v12', 'node': 3193, 'created_at': 'another installation'})
    result = plan_mapping(samples, proteins, tmp_path / 'plan', tmp_path / 'cache',
                          reference=marker, reference_pin=tmp_path / 'pin.json')
    assert result['reused_species'] == ['Alpha_plant']
    assert result['mapped_species'] == ['Beta_plant']


def test_same_release_reference_records_do_not_trigger_remapping(dataset_project):
    root = dataset_project
    imported(root)
    items = identities(root / 'input/metadata.tsv')[1]
    sources = [database_from_stages(root, root / 'resources/dataset_assets', members,
                                   name=name, reference_sha256=checksum * 64)
               for name, checksum, members in [('reference_a', 'a', items[:1]),
                                                ('reference_b', 'b', items[1:])]]
    config = root / 'config/build.yaml'
    cfg = yaml.safe_load(config.read_text())
    cfg['reuse_from'] = [str(p) for p in sources]
    config.write_text(yaml.safe_dump(cfg))
    report = dataset.plan(root, config)[-1]
    assert all(row[stage] == 'reuse' for row in report for stage in dataset.STAGES)
    assert all(row['mapping'] == 'reuse' for row in report)
    build = dataset.prepare(root, 'same_release', config)
    cached = load_tables(build / 'work/cache/odb/v12_3193/reused/snapshot.json')
    assert cached['reference_sha256s'] == ['a' * 64, 'b' * 64]


def test_known_snapshot_supersedes_unknown_cache_on_retry(tmp_path):
    proteins = tmp_path / 'proteins'
    proteins.mkdir()
    (proteins / 'Alpha_plant_protein.fa').write_text('>Alpha_plant_g1\nMK\n')
    snapshot(tmp_path / 'cache/v12_3193/aaa_unknown', proteins, ['Alpha_plant'])
    known = tmp_path / 'cache/v12_3193/zzz_known'
    snapshot(known, proteins, ['Alpha_plant'])
    data = json.loads((known / 'snapshot.json').read_text())
    data['reference_sha256'] = 'a' * 64
    write_json(known / 'snapshot.json', data)
    samples = tmp_path / 'samples.tsv'
    write_tsv(samples, ['species', 'odb_species'], [{'species': 'Alpha_plant', 'odb_species': 'Alpha_plant'}])
    marker = tmp_path / 'reference.json'
    write_json(marker, {'version': 'v12', 'node': 3193, 'created_at': 'another installation'})
    result = plan_mapping(samples, proteins, tmp_path / 'plan', tmp_path / 'cache',
                          reference=marker, reference_pin=tmp_path / 'pin.json')
    assert result['mapped_species'] == []
    assert result['sources'][0]['root'] == str(known)
    assert json.loads((tmp_path / 'pin.json').read_text())['reference_sha256'] is None


def test_reference_namespaces_are_distinct():
    assert reference_directory({'ncbi_tax_id': 3193}) == 'resources/orthodb/v12_3193'
    assert reference_directory({'version': 'v11', 'ncbi_tax_id': 3193}) == 'resources/orthodb/v11_3193'
