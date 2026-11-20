"""Compressed outputs retain expression values, alignment records, and export checks."""
import gzip
import json
import shutil

from alignment_batches import batch_targets, plan, run_batch
from analysis_validation import validate
from common import file_record, read_tsv, text_open, write_json, write_tsv
from filter_species import export, fasta_records
from og_kegg import aggregate
from test_alignment_batches import fixture
from test_filter_species import snapshot, SPECIES
from test_og_kegg import og_inputs, og_annotations, og_expression


def test_compressed_ko_tables_preserve_values_and_record_actual_files(og_inputs, og_annotations, og_expression, tmp_path):
    plain, compressed = tmp_path / 'plain', tmp_path / 'gzip'
    aggregate(og_inputs['samples'], og_annotations, og_expression, plain)
    aggregate(og_inputs['samples'], og_annotations, og_expression, compressed, compressed=True)
    for name in ('ko_tpm_sum.tsv', 'ko_tpm_sum_wide.tsv', 'ko_support.tsv'):
        with text_open(compressed / (name + '.gz')) as handle:
            assert handle.read() == (plain / name).read_text()
        assert not (compressed / name).exists()
    report = json.loads((compressed / 'provenance.json').read_text())
    assert {entry['path'] for entry in report['results']} == {
        str((compressed / name).resolve()) for name in
        ('ko_tpm_sum.tsv.gz', 'ko_tpm_sum_wide.tsv.gz', 'ko_support.tsv.gz', 'mapping_qc.tsv')}
    for use_gzip in (False, True):
        aggregate(og_inputs['samples'], og_annotations, og_expression, compressed, compressed=use_gzip)
        report = json.loads((compressed / 'provenance.json').read_text())
        assert {entry['path'] for entry in report['results']} == {
            str((compressed / (name + ('.gz' if use_gzip and name != 'mapping_qc.tsv' else ''))).resolve())
            for name in ('ko_tpm_sum.tsv', 'ko_tpm_sum_wide.tsv', 'ko_support.tsv', 'mapping_qc.tsv')}
        assert not list(compressed.glob('ko_*.tsv' if use_gzip else 'ko_*.tsv.gz'))


def test_gzip_batches_reuse_plain_cache_and_validate_compressed_members(tmp_path):
    inputs, plans, out, reports = fixture(tmp_path, count=2)
    cache = tmp_path / 'cache'
    plan(inputs, plans)
    batch, receipt = plans / 'batch_000000.json', out / 'batches/batch_000000.json'
    run_batch(batch, inputs, out, reports, receipt)
    cache.mkdir()
    for og in ('OG0', 'OG1'):
        (out / (og + '.faa')).replace(cache / (og + '.faa'))
        shutil.copyfile(reports / (og + '.json'), cache / (og + '.json'))
    run_batch(batch, inputs, out, reports, receipt, cache=cache, compressed=True)
    assert json.loads(receipt.read_text())['reused'] == 2
    assert len(batch_targets(plans, out, reports)) == 1 and receipt.exists()
    for og in ('OG0', 'OG1'):
        output = out / (og + '.faa.gz')
        assert list(fasta_records(output)) == list(fasta_records(inputs / (og + '.faa')))
        assert json.loads((reports / (og + '.json')).read_text())['alignment'] == file_record(output)
        assert not (out / (og + '.faa')).exists()
    run_batch(batch, inputs, out, reports, receipt, compressed=False)
    assert json.loads(receipt.read_text())['reused'] == 2
    assert not list(out.glob('*.faa.gz'))
    run_batch(batch, inputs, out, reports, receipt, compressed=True)
    assert json.loads(receipt.read_text())['reused'] == 2
    assert not list(out.glob('*.faa'))


def test_compressed_analysis_exports_and_verifies_without_changing_residues(snapshot):
    source, traits = snapshot
    metadata = source / 'metadata/metadata_high_busco.tsv'
    rows = [{**row, 'family': 'Plantaceae'} for row in read_tsv(metadata)]
    write_tsv(metadata, list(rows[0]), rows)
    paths = list((source / 'orthogroups/expression').glob('tpm*.tsv'))
    paths += list((source / 'kegg').glob('ko_tpm_sum*.tsv'))
    report_path = source / 'orthogroups/alignments/provenance.json'
    report = json.loads(report_path.read_text())
    for entry in report['alignments']:
        paths.append(source / 'orthogroups/alignments' / (entry['orthogroup'] + '.faa'))
    for path in paths:
        compressed = path.with_name(path.name + '.gz')
        compressed.write_bytes(gzip.compress(path.read_bytes(), mtime=0))
        path.unlink()
    for entry in report['alignments']:
        entry['alignment'] = file_record(source / 'orthogroups/alignments' / (entry['orthogroup'] + '.faa.gz'))
    write_json(report_path, report)
    filtered = source / 'filtered'
    export(source, ['Plant_A'], filtered, traits)
    validate(source, exclusions=['Plant_A'])
    rows = read_tsv(filtered / 'orthogroups/expression/tpm.tsv.gz')
    assert {r['species'] for r in rows} == set(SPECIES[1:])
    assert {r['tpm'] for r in rows} == {'1000000.000000'}
    assert (filtered / 'kegg/ko_tpm_sum.tsv.gz').is_file()
    alignment = filtered / 'orthogroups/alignments/OGshared.faa.gz'
    assert all(sequence == '--C--D' for _, sequence in fasta_records(alignment))
