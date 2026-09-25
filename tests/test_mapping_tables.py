"""Mapping reuse is species-local, bounded in memory, and independent of SQLite."""
import gzip
import json
from pathlib import Path

import pytest

from common import write_tsv
from incremental_odb import import_snapshot, plan
from mapping_tables import collect, load_tables, read_species, write_tables
from test_incremental_odb import snapshot


@pytest.fixture
def inputs(tmp_path):
    proteins = tmp_path/'proteins'; proteins.mkdir()
    names = ['Alpha_plant','Beta_plant','Removed_plant']
    for name in names:
        (proteins/f'{name}_protein.fa').write_text(f'>{name}_g1\nMK\n>{name}_g2\nMP\n>{name}_g3\nMM\n')
    old = tmp_path/'old'; snapshot(old, proteins, names)
    samples = tmp_path/'samples.tsv'
    write_tsv(samples, ['species','odb_species'], [{'species':s,'odb_species':s} for s in names])
    return proteins, old, samples, names


def test_tables_reuse_on_removal_and_readdition_without_parsing_inputs(inputs, tmp_path, monkeypatch):
    import mapping_tables
    proteins, old, samples, names = inputs
    cache = tmp_path/'cache'
    out = collect(samples,'unused','unused',proteins,tmp_path/'first',cache,existing=old)
    first = load_tables(out)
    assert read_species(out,names[0])[0][names[0]+'_g3'] == []
    def forbidden(*args, **kwargs): raise AssertionError('unchanged inputs must not be reparsed')
    monkeypatch.setattr(mapping_tables,'fasta_ids',forbidden)
    monkeypatch.setattr(mapping_tables,'annotation_pairs',forbidden)
    for label, selected in [('removed',names[:1]), ('restored',names)]:
        write_tsv(samples,['species','odb_species'],[{'species':s,'odb_species':s} for s in selected])
        result=collect(samples,'unused','unused',proteins,tmp_path/label,cache,existing=old)
        current=load_tables(result)
        assert set(current['tables']) == set(selected)
        assert {s:e['table']['sha256'] for s,e in current['tables'].items()} == {s:first['tables'][s]['table']['sha256'] for s in selected}
    assert not list(tmp_path.rglob('*.sqlite'))


def test_corrupt_species_cache_is_rejected(inputs, tmp_path):
    proteins,old,samples,names=inputs
    collect(samples,'unused','unused',proteins,tmp_path/'first',tmp_path/'cache',existing=old)
    cached=next((tmp_path/'cache/.tables').glob('*/table.tsv.gz'))
    value=bytearray(cached.read_bytes());value[-1]^=1;cached.write_bytes(value)
    with pytest.raises(ValueError,match='registered file changed'):
        collect(samples,'unused','unused',proteins,tmp_path/'second',tmp_path/'cache',existing=old)
    assert not (tmp_path/'second/snapshot.json').exists()


def test_changed_protein_maps_only_changed_species(inputs, tmp_path):
    proteins,old,samples,names=inputs
    cache=tmp_path/'cache';import_snapshot(old,cache)
    path=proteins/f'{names[0]}_protein.fa';path.write_text(path.read_text().replace('MK','ML'))
    result=plan(samples,proteins,tmp_path/'plan',cache)
    assert result['mapped_species'] == [names[0]]
    assert result['reused_species'] == names[1:]
