"""Expression reuse is bound to the sample, abundance and mapping."""
import json
from pathlib import Path

import pytest

import aggregate_tpm
from common import read_tsv, write_tsv
from mapping_fixtures import make_mapping, edit_mapping
from prepare_metadata import prepare


@pytest.mark.parametrize('change', ['abundance', 'mapping', 'ambiguous_mapping'])
def test_expression_cache_reuse_and_invalidation(tiny_inputs, tmp_path, monkeypatch, change):
    meta = tmp_path / 'metadata'
    prepare(**tiny_inputs, outdir=meta)
    mapping = tmp_path / 'mapping/snapshot.json'
    genes = [(f'Alpha_plant_g{i}', 'Alpha_plant') for i in (1, 2, 3)]
    make_mapping(mapping, genes, [('Alpha_plant_g1', 'OG1'), ('Alpha_plant_g2', 'OG2')])
    cache = tmp_path / 'cache'
    def run(name):
        out = tmp_path / name
        aggregate_tpm.aggregate(meta / 'samples.tsv', 'A1', mapping, out / 'tpm.tsv', out / 'qc.json', cache_dir=cache)
        return out
    first = run('first')
    implementation = aggregate_tpm._aggregate
    def unexpected(*args, **kwargs):
        pytest.fail('unchanged sample expression was recomputed')
    monkeypatch.setattr(aggregate_tpm, '_aggregate', unexpected)
    second = run('second')
    assert (first / 'tpm.tsv').samefile(second / 'tpm.tsv')
    monkeypatch.setattr(aggregate_tpm, '_aggregate', implementation)
    if change == 'abundance':
        abundance = Path(next(r for r in read_tsv(meta / 'samples.tsv') if r['run'] == 'A1')['abundance'])
        rows = read_tsv(abundance); rows[0]['tpm'] = '80'
        write_tsv(abundance, list(rows[0]), rows)
    elif change == 'mapping':
        edit_mapping(mapping, pairs=lambda rows: [(g, 'OG4' if og == 'OG1' else og) for g, og in rows])
    else:
        edit_mapping(mapping, pairs=lambda rows: rows + [('Alpha_plant_g1', 'OG3')])
        with pytest.raises(ValueError, match='A1: 1 genes map to multiple OGs; examples.*Alpha_plant_g1'):
            run('third')
        assert not (tmp_path / 'third/tpm.tsv').exists()
        assert not (tmp_path / 'third/qc.json').exists()
        assert len(list(cache.glob('*/receipt.json'))) == 1
        return
    third = run('third')
    assert not (first / 'tpm.tsv').samefile(third / 'tpm.tsv')
    assert (first / 'tpm.tsv').read_bytes() != (third / 'tpm.tsv').read_bytes()
    assert len(list(cache.glob('*/receipt.json'))) == 2
    assert json.loads((third / 'qc.json').read_text())['multimap'] == 'error'


def test_concurrent_builds_share_one_expression_calculation(tiny_inputs, tmp_path, monkeypatch):
    from concurrent.futures import ThreadPoolExecutor
    from threading import Barrier
    import time
    meta = tmp_path / 'metadata'
    prepare(**tiny_inputs, outdir=meta)
    mapping = tmp_path / 'mapping/snapshot.json'
    make_mapping(mapping, [('Alpha_plant_g1', 'Alpha_plant')], [('Alpha_plant_g1', 'OG1')])
    implementation = aggregate_tpm._aggregate
    calls = []
    def calculate(*args, **kwargs):
        calls.append(True)
        time.sleep(0.1)
        return implementation(*args, **kwargs)
    monkeypatch.setattr(aggregate_tpm, '_aggregate', calculate)
    barrier = Barrier(2)
    def run(name):
        barrier.wait()
        output = tmp_path / name / 'tpm.tsv'
        aggregate_tpm.aggregate(meta / 'samples.tsv', 'A1', mapping, output, output.with_suffix('.json'),
                                cache_dir=tmp_path / 'cache')
        return output
    with ThreadPoolExecutor(max_workers=2) as executor:
        first, second = list(executor.map(run, ['first', 'second']))
    assert calls == [True]
    assert first.samefile(second)
