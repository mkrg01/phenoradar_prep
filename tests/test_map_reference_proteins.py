import gzip
import json
import pytest

from common import file_record
import map_reference_proteins as reference_mapping
from map_reference_proteins import map_reference_proteins


def test_reference_mapping_retains_locus_ids(tmp_path, frozen_reference, monkeypatch):
    proteins = tmp_path / 'proteins.fa.gz'
    with gzip.open(proteins, 'wt') as handle:
        handle.write('>Sobic.003G036200 description\nMKK\n>Sobic.010G160700\nMKT\n>Sobic.001G000100\nMKI\n')
    annotation = ('#query\tODB_OG\nSobic.010G160700\t2at3193\nSobic.010G160700\t1at3193\n'
                  'Sobic.010G160700\t2at3193\nSobic.001G000100\t3at3193\n')
    reference = frozen_reference / 'reference.json'

    def run_odb(manifest, ref, output, work, label, **options):
        assert manifest.read_text() == f'{proteins.resolve()}\n'
        assert ref == reference.resolve()
        assert options['version'] == 'v12' and options['node'] == 3193
        assert options['gzip_results']
        output.mkdir(exist_ok=True)
        with gzip.open(output / f'{label}.og.annotations.gz', 'wt') as handle:
            handle.write(annotation)
        (output / 'provenance.json').write_text('{"fingerprint": "fake-mapping"}\n')
        return {'fingerprint': 'fake-mapping'}

    monkeypatch.setattr(reference_mapping, 'run_odb', run_odb)
    output = tmp_path / 'output'
    provenance = map_reference_proteins(proteins, reference, output)
    table = output / 'gene_to_orthogroup.tsv.gz'
    with gzip.open(table, 'rt') as handle:
        assert handle.read() == ('gene_id\torthogroup\nSobic.001G000100\t3at3193\nSobic.003G036200\t\n'
                                 'Sobic.010G160700\t1at3193\nSobic.010G160700\t2at3193\n')
    assert provenance == json.loads((output / 'mapping_provenance.json').read_text())
    assert provenance['table'] == file_record(table)
    assert provenance['orthodb']['version'] == 'v12'
    assert provenance['orthodb']['node'] == 3193
    assert provenance['qc'] == {'protein_genes': 3, 'mapped_genes': 2, 'unmapped_genes': 1,
                                'ambiguous_genes': 1, 'unique_gene_og_pairs': 3, 'duplicate_pairs_removed': 1}
    assert map_reference_proteins(proteins, reference, output) == provenance

    annotation = '#query\tODB_OG\nunknown\t3at3193\n'
    with pytest.raises(ValueError, match='query is absent'):
        map_reference_proteins(proteins, reference, tmp_path / 'bad-output')
    assert not (tmp_path / 'bad-output' / 'gene_to_orthogroup.tsv.gz').exists()

    with gzip.open(proteins, 'at') as handle:
        handle.write('>new_gene\nMKT\n')
    with pytest.raises(ValueError, match='mapping inputs changed'):
        map_reference_proteins(proteins, reference, output)
    assert file_record(table) == provenance['table']
