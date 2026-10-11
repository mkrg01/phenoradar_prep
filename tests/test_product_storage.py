"""Compression preserves scientific content, reuse, and stored-byte integrity."""
import gzip
import json
import shutil
from pathlib import Path

import pytest

from common import content_sha256, read_tsv, write_tsv
from dataset_assets import record, verify, verify_many
from product_storage import product_record, stage_product
from verification_cache import verification_cache
from test_datasets import dataset_project, imported


def test_compression_is_lossless_deterministic_and_shared(tmp_path):
    source = tmp_path / 'abundance.tsv'
    data = b'target_id\tlength\teff_length\test_counts\ttpm\r\nSpecies_g1\t9\t1.0000\t0\t0.000000\r\n' * 100
    source.write_bytes(data)
    original = record(source)
    first = stage_product(original, tmp_path / 'first.tsv.gz')
    second = stage_product(original, tmp_path / 'second.tsv.gz')
    assert Path(first['path']).read_bytes() == Path(second['path']).read_bytes()
    assert gzip.decompress(Path(first['path']).read_bytes()) == data
    assert first['bytes'] < original['bytes']
    assert content_sha256(first) == original['sha256']
    linked = stage_product(first, tmp_path / 'linked.tsv.gz')
    assert Path(linked['path']).samefile(first['path'])
    plain = stage_product(linked, tmp_path / 'decoded.tsv')
    assert Path(plain['path']).read_bytes() == data
    assert source.read_bytes() == data


def test_gzip_tsv_retains_numeric_text_zeros_and_columns(tmp_path):
    rows = [{'target_id': 'Species_g1', 'tpm': '0.000000', 'est_counts': '0'},
            {'target_id': 'Species_g2', 'tpm': '1.2345678901234567', 'est_counts': '3.5000'}]
    path = tmp_path / 'abundance.tsv.gz'
    write_tsv(path, list(rows[0]), rows)
    assert read_tsv(path) == rows
    first = path.read_bytes()
    write_tsv(path, list(rows[0]), rows)
    assert path.read_bytes() == first


def test_moved_compressed_product_verifies_stored_bytes(tmp_path):
    source = tmp_path / 'protein.fa'
    source.write_text('>Species_g1\nMK*\n' * 100)
    compressed = stage_product(record(source), tmp_path / 'protein.fa.gz')
    moved = tmp_path / 'moved.fa.gz'
    moved.write_bytes(Path(compressed['path']).read_bytes())
    entry = dict(compressed, path=str(moved))
    with verification_cache(tmp_path / 'checks.sqlite'):
        verify(entry)
        assert verify_many([entry]) == [None]
        moved.write_bytes(moved.read_bytes() + b'changed')
        with pytest.raises(ValueError, match='registered file changed'):
            verify(entry)
        assert isinstance(verify_many([entry])[0], ValueError)


def test_conversion_rejects_changed_source_and_conflicting_destination(tmp_path):
    source = tmp_path / 'protein.fa'
    source.write_text('>Species_g1\nMK*\n')
    original = record(source)
    source.write_text('>Species_g1\nML*\n')
    with pytest.raises(ValueError, match='(?:registered file|product content) changed'):
        stage_product(original, tmp_path / 'changed.fa.gz')
    assert not (tmp_path / 'changed.fa.gz').exists()
    destination = tmp_path / 'existing.fa.gz'
    with gzip.open(destination, 'wt') as handle:
        handle.write('different protein\n')
    before = destination.read_bytes()
    with pytest.raises(ValueError, match='refusing to replace'):
        stage_product(record(source), destination)
    assert destination.read_bytes() == before


def legacy_copy(source, target):
    """Reconstruct the historical schema-5 plain storage layout for compatibility."""
    from common import write_json
    from dataset_assets import digest
    shutil.copytree(source, target)
    manifest = target / 'manifest.json'
    data = json.loads(manifest.read_text())
    inventory = {entry['path']: entry for entry in data['files']}

    def replace(path):
        relative = str(path.relative_to(target))
        entry = dict(record(path), path=relative)
        inventory[relative] = entry
        return entry

    for product in data['products'].values():
        for key in ('abundance', 'protein', 'expression'):
            compressed = target / product[key]['path']
            plain = compressed.with_suffix('')
            plain.write_bytes(gzip.decompress(compressed.read_bytes()))
            inventory.pop(product[key]['path'])
            compressed.unlink()
            product[key] = replace(plain)
        translation = target / product['translation']['path']
        info = json.loads(translation.read_text())
        info['protein'] = product['protein']
        write_json(translation, info)
        product['translation'] = replace(translation)
        qc = target / product['expression_qc']['path']
        info = json.loads(qc.read_text())
        info.update(abundance=product['abundance'], expression=product['expression'])
        write_json(qc, info)
        product['expression_qc'] = replace(qc)
    data['files'] = list(inventory.values())
    data['sha256'] = digest({key: value for key, value in data.items() if key != 'sha256'})
    write_json(manifest, data)
    return target


def test_legacy_and_compressed_databases_share_reuse_and_imports(dataset_project):
    import dataset
    import yaml
    from build_products import import_expression, import_protein, load_complete
    from common import write_json
    from mapping_tables import load_tables, protein_record

    root = dataset_project
    compressed = imported(root)
    legacy = legacy_copy(compressed, root / 'resources/legacy/database')
    original = {p: p.read_bytes() for p in legacy.rglob('*') if p.is_file()}
    old = load_complete(legacy)
    modern = load_complete(compressed)
    for name in old['products']:
        for key in ('abundance', 'protein', 'expression'):
            assert old['products'][name][key]['sha256'] == content_sha256(modern['products'][name][key])
    config = root / 'config/build.yaml'
    cfg = yaml.safe_load(config.read_text())
    cfg['reuse_from'] = [str(legacy), str(compressed)]
    config.write_text(yaml.safe_dump(cfg))
    report = dataset.plan(root, config)[-1]
    assert all(row['mapping'] == 'reuse' for row in report)
    build = dataset.prepare(root, 'mixed_storage', config)
    inputs = dataset.materialize(build)
    assert len(list(inputs.glob('quant/*/*/*_abundance.tsv.gz'))) == len(old['products'])
    protein = build / 'test/Alpha_plant_A1_protein.fa.gz'
    provenance = protein.with_suffix('').with_suffix('.json')
    import_protein(legacy, 'Alpha_plant_A1', protein, provenance)
    expected = load_tables(old['mapping']['path'])['tables']['Alpha_plant_A1']['protein_sha256']
    assert content_sha256(protein_record(protein)) == expected
    product = old['products']['Alpha_plant_A1']
    samples = build / 'test/samples.tsv'
    write_tsv(samples, ['species', 'run', 'abundance'],
              [dict(species='Alpha_plant_A1', run='A1', abundance=product['abundance']['path'])])
    expression, qc = build / 'test/A1.tsv.gz', build / 'test/A1.qc.json'
    import_expression(legacy, samples, 'A1', expression, qc)
    assert read_tsv(expression) == read_tsv(product['expression']['path'])
    assert json.loads(qc.read_text())['expression']['sha256'] == record(expression)['sha256']
    assert all(p.read_bytes() == data for p, data in original.items())
