import json

import pytest

from common import write_json
from dataset_assets import record
from upgrade_analysis import verified_imports


def test_upgrade_only_reuses_source_equivalent_imports(tmp_path):
    source = tmp_path / 'source'; source.mkdir()
    protein = source / 'protein.fa.gz'; protein.write_bytes(b'protein product')
    expression = source / 'expression.tsv.gz'; expression.write_bytes(b'expression product')
    translation = source / 'translation.json'
    qc = source / 'qc.json'
    write_json(translation, {'protein': record(protein), 'translation_table': 1})
    write_json(qc, {'expression': record(expression), 'run': 'R1'})
    product = {'odb_species': 'Plant_R1', 'row': {'run': 'R1'},
               'protein': record(protein), 'expression': record(expression),
               'translation': record(translation), 'expression_qc': record(qc)}
    out = tmp_path / 'analysis'
    paths = [out / 'proteins/Plant_R1_protein.fa.gz', out / 'proteins/Plant_R1_protein.json',
             out / 'orthogroups/expression/runs/R1.tsv.gz', out / 'orthogroups/expression/runs/R1.qc.json']
    for target, original in zip(paths, [protein, translation, expression, qc]):
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(original.read_bytes())
    manifest = {'root': str(tmp_path), 'pipeline': {'output_root': 'analysis'}}
    completed = {'products': {'Plant_R1': product}}
    assert len(verified_imports(out, manifest, completed)) == 4
    data = json.loads(paths[1].read_text()); data['translation_table'] = 2
    write_json(paths[1], data)
    with pytest.raises(ValueError, match='changed imported provenance'):
        verified_imports(out, manifest, completed)
