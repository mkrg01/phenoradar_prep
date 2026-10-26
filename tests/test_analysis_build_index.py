import json

import pytest

from analysis_build_index import BuildInputs, create_index
from common import write_json
from dataset_assets import record


def test_lazy_workers_do_not_load_large_database_and_shards_are_verified(tmp_path, monkeypatch):
    source = tmp_path / 'database.json'; source.write_text('{}')
    import portable_build
    monkeypatch.setattr(portable_build, 'completion_path', lambda p: source)
    product = {k: record(source) for k in ('cds', 'abundance', 'protein', 'translation', 'expression', 'expression_qc')}
    product.update(odb_species='Plant_R1', row={'run': 'R1', 'large_metadata': 'x' * 10000})
    data = {'products': {'Plant_R1': product}, 'mapping': record(source), 'tpm': {'multimap': 'error'}}
    index = create_index(source, tmp_path / 'index', data)
    import build_products
    monkeypatch.setattr(build_products, 'load_complete', lambda *a, **k: pytest.fail('large DB read on worker'))
    lazy = BuildInputs(source, index)
    assert lazy._data is None
    loaded = lazy.product('R1', 'run')
    assert loaded['row'] == {'run': 'R1'}
    assert index.stat().st_size < 2000
    shard = next((index.parent / 'products').glob('*.json'))
    shard.write_text(shard.read_text().replace('Plant_R1', 'Plant_R2'))
    with pytest.raises(ValueError, match='registered file changed'):
        BuildInputs(source, index).product('R1', 'run')
    source.write_text('changed database')
    with pytest.raises(ValueError, match='registered file changed'):
        BuildInputs(source, index).data
