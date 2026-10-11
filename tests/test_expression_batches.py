"""Batched jobs preserve run outputs, validation, and independent retry checkpoints."""
import builtins
import json
from pathlib import Path
import shutil

import pytest

import aggregate_tpm as expression
import mapping_tables
import sample_table
from common import write_tsv
from mapping_fixtures import make_mapping


@pytest.fixture
def inputs(tmp_path):
    rows, genes, pairs = [], [], []
    for species, run in [('Plant_A', 'A1'), ('Plant_B', 'B1'), ('Plant_C', 'C1')]:
        abundance = tmp_path / f'{run}.abundance.tsv'
        write_tsv(abundance, ['target_id', 'tpm'],
                  [{'target_id': f'{species}_g{i}', 'tpm': value}
                   for i, value in enumerate([20, 30, 50], 1)])
        rows.append(dict(species=species, odb_species=species, run=run, abundance=str(abundance)))
        genes.extend((f'{species}_g{i}', species) for i in range(1, 4))
        pairs.extend((f'{species}_g{i}', f'OG{i}') for i in (1, 2))
    samples = tmp_path / 'samples.tsv'
    write_tsv(samples, list(rows[0]), rows)
    mapping = make_mapping(tmp_path / 'mapping/snapshot.json', genes, pairs)
    return samples, mapping, rows


def comparable(value):
    if isinstance(value, dict):
        return {key: comparable(item) for key, item in value.items()
                if key not in {'created_at', 'path', 'stat'}}
    if isinstance(value, list):
        return [comparable(item) for item in value]
    return value


@pytest.mark.parametrize('compressed', [False, True])
def test_batch_matches_single_runs_and_reads_shared_indexes_once(inputs, tmp_path, monkeypatch, compressed):
    samples, mapping, rows = inputs
    suffix = '.tsv.gz' if compressed else '.tsv'
    for row in rows:
        run = row['run']
        expression.aggregate(samples, run, mapping, tmp_path / 'single' / (run + suffix),
                             tmp_path / 'single' / f'{run}.qc.json', cache_dir=tmp_path / 'single-cache')
    sample_table._read.cache_clear()
    mapping_tables._read_tables.cache_clear()
    reads = []
    open_file, read_text = builtins.open, Path.read_text
    def counted_open(path, *args, **kwargs):
        if Path(path) == samples:
            reads.append('samples')
        return open_file(path, *args, **kwargs)
    def counted_text(path, *args, **kwargs):
        if path == mapping:
            reads.append('mapping')
        return read_text(path, *args, **kwargs)
    monkeypatch.setattr(builtins, 'open', counted_open)
    monkeypatch.setattr(Path, 'read_text', counted_text)
    batch = tmp_path / 'batch'
    expression.aggregate_many(samples, [r['run'] for r in rows], mapping, batch, compressed=compressed)
    assert reads == ['samples', 'mapping']
    for row in rows:
        run = row['run']
        assert (batch / (run + suffix)).read_bytes() == (tmp_path / 'single' / (run + suffix)).read_bytes()
        assert comparable(json.loads((batch / f'{run}.qc.json').read_text())) == comparable(
            json.loads((tmp_path / 'single' / f'{run}.qc.json').read_text()))
        expression.publish_run(batch, run, tmp_path / 'published' / (run + suffix),
                               tmp_path / 'published' / f'{run}.qc.json')
        assert (batch / (run + suffix)).samefile(tmp_path / 'published' / (run + suffix))


def test_failed_run_does_not_stop_others_and_retry_reuses_successes(inputs, tmp_path, monkeypatch):
    samples, mapping, rows = inputs
    failed = Path(rows[1]['abundance'])
    valid = failed.read_bytes()
    failed.write_text('target_id\ttpm\nPlant_B_g1\t-1\n')
    calls = []
    calculate = expression._aggregate
    def counted(samples, run, *args, **kwargs):
        calls.append(run)
        return calculate(samples, run, *args, **kwargs)
    monkeypatch.setattr(expression, '_aggregate', counted)
    batch = tmp_path / 'batch'
    with pytest.raises(ValueError, match='expression batch failed:'):
        expression.aggregate_many(samples, ['B1', 'A1', 'C1'], mapping, batch)
    assert calls == ['B1', 'A1', 'C1']
    assert not (batch / 'receipt.json').exists()
    assert len(list((tmp_path / '.expression-cache').glob('*/receipt.json'))) == 2
    shutil.rmtree(batch)  # Snakemake removes the failed directory output.
    failed.write_bytes(valid)
    calls.clear()
    expression.aggregate_many(samples, ['B1', 'A1', 'C1'], mapping, batch)
    assert calls == ['B1']
    assert json.loads((batch / 'receipt.json').read_text())['runs'] == ['B1', 'A1', 'C1']


def test_publication_rejects_corrupt_batch_before_exposing_outputs(inputs, tmp_path):
    samples, mapping, _ = inputs
    batch = tmp_path / 'batch'
    expression.aggregate_many(samples, ['A1'], mapping, batch)
    (batch / 'A1.qc.json').write_text('damaged\n')
    with pytest.raises(ValueError, match='registered file changed'):
        expression.publish_run(batch, 'A1', tmp_path / 'published.tsv', tmp_path / 'published.qc.json')
    assert not (tmp_path / 'published.tsv').exists()


@pytest.mark.parametrize('runs', [[], ['A1', 'A1'], ['../A1'], ['missing']])
def test_invalid_batch_selection_is_rejected(inputs, tmp_path, runs):
    samples, mapping, _ = inputs
    with pytest.raises(ValueError):
        expression.aggregate_many(samples, runs, mapping, tmp_path / 'batch')
    assert not (tmp_path / 'batch').exists()


def test_batch_plan_is_bounded_and_tracks_changed_selection(inputs):
    samples, _, rows = inputs
    many = [dict(rows[0], run=f'R{i:03d}') for i in range(101)]
    write_tsv(samples, list(many[0]), reversed(many))
    batches, by_run = sample_table.expression_batches(samples)
    assert [len(batch) for batch in batches.values()] == [50, 50, 1]
    assert by_run['R100'] == 'batch_00002'
    write_tsv(samples, list(many[0]), many[:1])
    batches, by_run = sample_table.expression_batches(samples)
    assert list(batches) == ['batch_00000']
    assert by_run == {'R000': 'batch_00000'}

