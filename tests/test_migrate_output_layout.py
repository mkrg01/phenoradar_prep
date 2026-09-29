"""Offline migration preserves expression values and rejects mismatched provenance."""
import gzip
import json
import os
import shutil
from pathlib import Path

import pytest

from build_products import load_complete
from common import write_json
from dataset_assets import digest, record
from migrate_output_layout import migrate


@pytest.fixture
def legacy_database(tmp_path):
    source = tmp_path / 'old/products'; source.mkdir(parents=True)
    run = tmp_path / 'run001'
    receipts = tmp_path / 'receipts'; receipts.mkdir()
    species, old, accession = 'Alpha_plant_A1', 'Alpha_plant', 'A1'
    inventory = {}
    def add(relative, text):
        path = source / relative; path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
        entry = dict(record(path), path=relative)
        inventory[relative] = entry
        return entry
    product = {'row': {'scientific_name': 'Alpha plant', 'run': accession, 'taxid': '1',
                       'species_id': old, 'analysis_sample_id': species}, 'odb_species': species}
    for key, relative in [('cds', 'cds/sample.fa.gz'), ('busco', 'busco/full/sample.tsv'),
                          ('abundance', 'quant/sample.tsv'), ('protein', 'proteins/sample.fa'),
                          ('translation', 'proteins/sample.json')]:
        product[key] = add(relative, key + '\n')
    for relative in ('metadata.tsv', 'busco/summary.tsv', 'excluded_accessions.tsv', 'excluded_runs.tsv'):
        add(relative, 'header\n')
    table = source / f'odb/species/{species}.tsv.gz'; table.parent.mkdir(parents=True)
    with gzip.open(table, 'wt') as f:
        f.write(f'gene_id\torthogroup\n{species}@g1\tOG1\n{species}@g2\tOG2\n')
    inventory[str(table.relative_to(source))] = dict(record(table), path=str(table.relative_to(source)))
    table_record = dict(record(table), path=f'species/{species}.tsv.gz')
    annotation_sha = 'a' * 64
    entry = {'odb_species': species, 'protein_sha256': product['protein']['sha256'],
             'table': table_record, 'qc': {'protein_genes': 2}, 'annotation_sha256': annotation_sha}
    mapping = {'schema_version': 2, 'kind': 'odb_tables', 'version': 'v12', 'node': 3193,
               'proteins': [{'species': species, 'sha256': product['protein']['sha256']}],
               'tables': {species: entry}}
    write_json(source / 'odb/snapshot.json', mapping)
    mapping_record = dict(record(source / 'odb/snapshot.json'), path='odb/snapshot.json')
    inventory[mapping_record['path']] = mapping_record
    data = {'schema_version': 4, 'kind': 'completed_build', 'build_id': 'old', 'created_at': 'historical',
            'fields': list(product['row']), 'products': {species: product}, 'files': list(inventory.values()),
            'mapping': mapping_record, 'translation': {'table': 1}, 'lineage': 'embryophyta_odb12',
            'odb': {'node': 3193, 'version': 'v12'}, 'excluded_runs': []}
    data['sha256'] = digest(data); write_json(source / 'manifest.json', data)
    old_abundance, old_protein, old_mapping = 'b' * 64, 'c' * 64, 'd' * 64
    snapshot = tmp_path / 'legacy_snapshot.json'
    write_json(snapshot, {'schema_version': 1, 'node': 3193, 'version': 'v12',
                         'annotations': {'sha256': annotation_sha},
                         'proteins': [{'species': old, 'sha256': old_protein}]})
    write_json(run / 'orthogroups/mapping/merge_qc.json', {'sources': [record(snapshot)]})
    write_json(run / 'run.json', {'config': {'translation': data['translation'], 'odb': data['odb'],
               'phylogeny': {'lineage': data['lineage']}, 'tpm': {'multimap': 'error'}}})
    qc_path = run / 'orthogroups/expression/runs/A1.qc.json'
    write_json(qc_path, {'species': old, 'run': accession, 'multimap': 'error', 'protein_genes': 2,
                        'orthogroups': 2, 'retained_tpm': 5.0, 'abundance': {'sha256': old_abundance}})
    qc_path.with_suffix('').with_suffix('.tsv').write_text(
        'species\trun\torthogroup\ttpm_sum\ttpm\n'
        'Alpha_plant\tA1\tOG1\t2.0000\t400000.0\nAlpha_plant\tA1\tOG2\t3.0\t600000.00\n')
    conversions = {}
    for key, old_sha, new in [('quant', old_abundance, product['abundance']),
                              ('protein', old_protein, product['protein']),
                              ('mapping', old_mapping, table_record)]:
        conversions[key] = {'old_id': old, 'new_id': species,
                            'source': {'sha256': old_sha}, 'output': new}
    write_json(receipts / f'{species}.json', {'sample_id': species, 'old_id': old, 'run': accession,
               'conversions': conversions, 'mapping': {'annotation_sha256': annotation_sha}})
    return source, run, receipts


def test_migration_is_portable_and_preserves_numeric_text(legacy_database, tmp_path):
    source, run, receipts = legacy_database
    before = (source / 'manifest.json').read_bytes()
    old_expression = (run / 'orthogroups/expression/runs/A1.tsv').read_bytes()
    target = tmp_path / 'results/latest'
    migrate(source, target, run, receipts, workers=2)
    data = load_complete(target)
    assert data['schema_version'] == 5
    assert data['tpm'] == {'multimap': 'error'}
    expression = target / 'database/expression/runs/A1.tsv'
    assert expression.read_bytes() == old_expression.replace(b'Alpha_plant\t', b'Alpha_plant_A1\t')
    assert (source / 'manifest.json').read_bytes() == before
    assert (run / 'orthogroups/expression/runs/A1.tsv').read_bytes() == old_expression
    assert os.stat(source / 'proteins/sample.fa').st_ino == os.stat(target / 'database/proteins/sample.fa').st_ino
    # The database alone remains usable after relocation and removal of its inputs.
    copied = tmp_path / 'elsewhere/database'; shutil.copytree(target / 'database', copied)
    source.parent.rename(tmp_path / 'offline-old')
    run.rename(tmp_path / 'offline-run')
    relocated = load_complete(copied)
    assert len(relocated['products']) == 1
    assert all(Path(e['path']).is_relative_to(copied) for e in relocated['files'])
    report = json.loads((target / 'logs/migration.json').read_text())
    assert report['scientific_recomputation'] is False
    with pytest.raises(ValueError, match='already exists'):
        migrate(copied, target, run, receipts)


@pytest.mark.parametrize('mismatch', ['abundance', 'mapping', 'normalization', 'duplicate'])
def test_rejects_mismatched_or_corrupt_expression(legacy_database, tmp_path, mismatch):
    source, run, receipts = legacy_database
    qc_path = run / 'orthogroups/expression/runs/A1.qc.json'
    if mismatch == 'abundance':
        qc = json.loads(qc_path.read_text()); qc['abundance']['sha256'] = 'bad'; write_json(qc_path, qc)
    elif mismatch == 'mapping':
        path = receipts / 'Alpha_plant_A1.json'
        receipt = json.loads(path.read_text()); receipt['mapping']['annotation_sha256'] = 'bad'; write_json(path, receipt)
    else:
        path = qc_path.with_suffix('').with_suffix('.tsv')
        text = path.read_text()
        text = text.replace('600000.00', '1') if mismatch == 'normalization' else text.replace('OG2', 'OG1')
        path.write_text(text)
    target = tmp_path / 'results/latest'
    with pytest.raises(ValueError):
        migrate(source, target, run, receipts)
    assert not target.exists()
    assert not list(target.parent.glob('.*.migrate-*')) if target.parent.exists() else True
    assert (source / 'manifest.json').exists()
