"""Input indexes must preserve metadata and never freeze mutable reuse decisions."""
import copy
import json
import os
from pathlib import Path

import pytest
import yaml

import dataset_assets
import portable_build
from build_products import load_complete
from common import write_json, write_tsv
from database_reuse import software_provenance
from dataset import plan
from dataset_assets import digest, identities, record
from test_datasets import dataset_project, imported
from verification_cache import verification_cache


def test_verified_product_json_index_is_reused_by_deferred_preparation(dataset_project, monkeypatch):
    import database_reuse
    root = dataset_project
    imported(root)
    config = root / 'config/build.yaml'
    cache = root / '.cache/verification.sqlite'
    with verification_cache(cache):
        expected = plan(root, config)[-1]
    with verification_cache(cache), monkeypatch.context() as cached:
        cached.setattr(database_reuse, 'product_json',
                       lambda *_, **__: pytest.fail('unchanged product JSON should use its index'))
        assert plan(root, config, verify_files=False)[-1] == expected


def test_deferred_product_json_index_detects_same_size_edit_with_restored_mtime(dataset_project):
    root = dataset_project
    source = imported(root)
    config = root / 'config/build.yaml'
    cache = root / '.cache/verification.sqlite'
    with verification_cache(cache):
        plan(root, config)
    product = next(iter(load_complete(source, verify_files=False)['products'].values()))
    qc = Path(product['expression_qc']['path'])
    previous = qc.stat()
    original = qc.read_bytes()
    changed = original.replace(b'"multimap": "error"', b'"multimap": "other"')
    assert changed != original and len(changed) == len(original)
    qc.write_bytes(changed)
    os.utime(qc, ns=(previous.st_atime_ns, previous.st_mtime_ns))
    with verification_cache(cache):
        report = plan(root, config, verify_files=False, repair=True)[-1]
    conflicts = [row for row in report if row['assembly'] == 'conflict']
    assert len(conflicts) == 1
    assert 'expression identity differs' in conflicts[0]['reason']
    with verification_cache(cache):
        report = plan(root, config, repair=True)[-1]
    assert sum(bool(row.get('rebuild')) for row in report) == 1


def save_manifest(path, data):
    data['sha256'] = digest({k: v for k, v in data.items() if k != 'sha256'})
    write_json(path, data)


def test_metadata_index_restores_columns_values_and_independent_rows(tmp_path, monkeypatch):
    metadata = tmp_path / 'metadata.tsv'
    fields = ['scientific_name', 'run', 'taxid', 'future', 'spaces']
    write_tsv(metadata, fields, [{'scientific_name': 'Alpha plant', 'run': 'A1',
                                'taxid': '42.0', 'future': '', 'spaces': '  '}])
    original = metadata.read_bytes()
    with verification_cache(tmp_path / 'cache.sqlite'):
        expected = identities(metadata, source=record(metadata))
    with verification_cache(tmp_path / 'cache.sqlite'), monkeypatch.context() as cached:
        cached.setattr(dataset_assets, '_read_identities',
                       lambda *_: pytest.fail('unchanged metadata should use its index'))
        fields, items = identities(metadata, source=record(metadata))
        assert (fields, items) == expected
        assert items[0]['row']['future'] == ''
        assert items[0]['row']['spaces'] == '  '
        assert items[0]['row']['taxid'] == '42'
        items[0]['row']['future'] = 'caller mutation'
        assert identities(metadata, source=record(metadata)) == expected
    assert metadata.read_bytes() == original


@pytest.mark.parametrize('change', ['fill_empty', 'add_column', 'drop_column', 'add_sample'])
def test_metadata_updates_rebuild_the_index(tmp_path, change):
    metadata = tmp_path / 'metadata.tsv'
    fields = ['scientific_name', 'run', 'taxid', 'future']
    rows = [{'scientific_name': 'Alpha plant', 'run': 'A1', 'taxid': '42', 'future': ''}]
    cache = tmp_path / 'cache.sqlite'
    write_tsv(metadata, fields, rows)
    with verification_cache(cache):
        identities(metadata, source=record(metadata))
    if change == 'fill_empty':
        rows[0]['future'] = 'new annotation'
    elif change == 'add_column':
        fields.append('new_column')
        rows[0]['new_column'] = 'new value'
    elif change == 'drop_column':
        fields.remove('future')
        rows[0].pop('future')
    else:
        rows.append(dict(rows[0], run='A2', future='new sample'))
    write_tsv(metadata, fields, rows)
    with verification_cache(cache):
        assert identities(metadata, source=record(metadata)) == dataset_assets._read_identities(metadata)


def test_metadata_index_rejects_same_size_edit_with_restored_mtime(tmp_path):
    metadata = tmp_path / 'metadata.tsv'
    fields = ['scientific_name', 'run', 'taxid']
    write_tsv(metadata, fields, [{'scientific_name': 'Alpha plant', 'run': 'A1', 'taxid': '42'}])
    cache = tmp_path / 'cache.sqlite'
    with verification_cache(cache):
        identities(metadata, source=record(metadata))
    before = metadata.stat()
    metadata.write_bytes(metadata.read_bytes().replace(b'42', b'xx'))
    os.utime(metadata, ns=(before.st_atime_ns, before.st_mtime_ns))
    with verification_cache(cache), pytest.raises(ValueError, match='positive taxid'):
        identities(metadata, source=record(metadata))


def test_manifest_index_preserves_missing_keys_and_refreshes_filled_columns(dataset_project, monkeypatch):
    root = dataset_project
    source = imported(root)
    manifest = source / 'manifest.json'
    data = json.loads(manifest.read_text())
    names = list(data['products'])
    data['fields'].extend(['future', 'spaces'])
    data['products'][names[0]]['row'].update(future='', spaces='  ')
    save_manifest(manifest, data)
    cache = root / '.cache/verification.sqlite'
    with verification_cache(cache):
        first = load_complete(source, verify_files=False)
    portable_build._MANIFESTS.clear()
    with verification_cache(cache), monkeypatch.context() as cached:
        read_text = Path.read_text
        def no_manifest_read(path, *args, **kwargs):
            if path == manifest:
                pytest.fail('unchanged manifest should use its validated index')
            return read_text(path, *args, **kwargs)
        cached.setattr(Path, 'read_text', no_manifest_read)
        cached.setattr(portable_build, 'digest', lambda *_: pytest.fail('manifest was already validated'))
        second = load_complete(source, verify_files=False)
        assert second == first
        assert second['products'][names[0]]['row']['future'] == ''
        assert 'future' not in second['products'][names[1]]['row']
        second['products'][names[0]]['row']['future'] = 'caller mutation'
        assert load_complete(source, verify_files=False) == first
    data['products'][names[0]]['row']['future'] = 'filled later'
    data['fields'].append('added_later')
    data['products'][names[0]]['row']['added_later'] = 'new column'
    save_manifest(manifest, data)
    with verification_cache(cache):
        updated = load_complete(source, verify_files=False)
        assert updated['products'][names[0]]['row'] == data['products'][names[0]]['row']
        assert updated['fields'] == data['fields']
    data['build_id'] = 'changed without updating the manifest checksum'
    write_json(manifest, data)
    with verification_cache(cache), pytest.raises(ValueError, match='build completion record changed'):
        load_complete(source, verify_files=False)


def test_cached_manifest_still_checks_symlink_targets(dataset_project):
    root = dataset_project
    source = imported(root)
    cache = root / '.cache/verification.sqlite'
    with verification_cache(cache):
        data = load_complete(source, verify_files=False)
    protein = Path(next(iter(data['products'].values()))['protein']['path'])
    outside = root / 'outside.fa'
    outside.write_bytes(protein.read_bytes())
    protein.unlink()
    protein.symlink_to(outside)
    with verification_cache(cache), pytest.raises(ValueError, match='escapes bundle'):
        load_complete(source, verify_files=False)


def test_cached_inputs_recheck_build_settings_and_exclusions(dataset_project):
    root = dataset_project
    imported(root)
    config = root / 'config/build.yaml'
    cfg = yaml.safe_load(config.read_text())
    cache = root / '.cache/verification.sqlite'
    with verification_cache(cache):
        assert {r['assembly'] for r in plan(root, config)[-1]} == {'reuse'}
    changed = copy.deepcopy(cfg)
    changed['translation']['table'] = 2
    config.write_text(yaml.safe_dump(changed))
    with verification_cache(cache):
        report = plan(root, config)[-1]
        assert {r['assembly'] for r in report} == {'conflict'}
        assert all('genetic code or BUSCO lineage differs' in r['reason'] for r in report)
    policy = root / 'config/new_exclusions.tsv'
    write_tsv(policy, ['accession', 'reason'], [{'accession': 'A1', 'reason': 'new exclusion'}])
    cfg['excluded_accessions'] = str(policy)
    config.write_text(yaml.safe_dump(cfg))
    with verification_cache(cache):
        report = plan(root, config)[-1]
        assert sum(r['assembly'] == 'reuse' for r in report) == 2
        assert sum(r['assembly'] == 'excluded' for r in report) == 1


def test_software_index_tracks_changed_pins_and_verifies_the_source(tmp_path, monkeypatch):
    source = tmp_path / 'build.json'
    build = {'config': {'genegalleon': {'version': '1.0', 'revision': 'abc', 'image_uri': 'old'}},
             'software_lock': {'source': 'original'}}
    write_json(source, build)
    entry = record(source)
    cache = tmp_path / 'cache.sqlite'
    with verification_cache(cache):
        expected = software_provenance(entry)
    with verification_cache(cache), monkeypatch.context() as cached:
        cached.setattr(Path, 'read_text', lambda *_: pytest.fail('historical metadata should not be reparsed'))
        assert software_provenance(entry) == expected
    build['config']['genegalleon']['image_uri'] = 'new'
    build['software_lock'] = {'source': 'updated'}
    write_json(source, build)
    with verification_cache(cache):
        with pytest.raises(ValueError, match='registered file changed'):
            software_provenance(entry)
        updated = software_provenance(record(source))
        assert updated['pins']['image_uri'] == 'new'
        assert updated['software_lock'] == {'source': 'updated'}
