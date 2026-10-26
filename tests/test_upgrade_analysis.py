import json
import copy
import subprocess

import yaml

import pytest

from common import write_json
from dataset_assets import record
from upgrade_analysis import verified_imports


@pytest.mark.parametrize('saved,requested,expected', [({}, None, (0, 0)), ({}, 0.5, (0, 0.5)),
    ({'min_sample_fraction': 0.75}, None, (0.75, 0.75)),
    ({'min_sample_fraction': 0.75}, 0, (0.75, 0))])
def test_occupancy_migration_requires_explicit_condition_change(saved, requested, expected):
    from upgrade_analysis import occupancy_condition
    assert occupancy_condition({'pipeline': {'alignment': saved}}, requested) == expected


@pytest.mark.parametrize('requested', [None, 0.5])
def test_upgrade_backs_up_and_refreezes_only_declared_condition(tmp_path, monkeypatch, requested):
    import upgrade_analysis as module
    from dataset import implementation
    from dataset_assets import verify
    root = tmp_path / 'project'; root.mkdir()
    (root / 'workflow/scripts').mkdir(parents=True)
    (root / 'workflow/scripts/align_orthogroups.py').write_text('unchanged scientific aligner\n')
    (root / 'workflow/Snakefile').write_text('workflow\n')
    (root / 'run_pipeline.sh').write_text('launcher\n')
    (root / 'VERSION').write_text('0.4.14\n')
    subprocess.run(['git', 'init', '-q', str(root)], check=True)
    def commit():
        subprocess.run(['git', 'add', '.'], cwd=root, check=True)
        subprocess.run(['git', '-c', 'user.name=Test', '-c', 'user.email=test@example.invalid',
                        'commit', '-qm', 'snapshot'], cwd=root, check=True)
        return subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=root, text=True).strip()
    monkeypatch.setattr(module, 'PREDECESSOR', commit())
    (root / 'VERSION').write_text('0.4.15\n')
    monkeypatch.setattr(module, 'BATCHED', commit())
    old_code = implementation(root)
    (root / 'VERSION').write_text('0.4.16\n')
    (root / 'workflow/scripts/alignment_selection.py').write_text('occupancy selection\n')
    commit()
    path = root / 'results/test'; (path / 'jobs').mkdir(parents=True)
    config = {'alignment': {'enabled': True}, 'unchanged': {'value': 17}}
    for filename in ('analysis.yaml', 'pipeline.yaml'):
        (path / filename).write_text(yaml.safe_dump(config))
    previous_config = (path / 'pipeline.yaml').read_bytes()
    output = path / 'tpm.tsv'; output.write_text('saved TPM\n')
    old_output = record(output)
    source = root / 'database.json'; source.write_text('{}')
    manifest = {'root': str(root), 'pipeline': copy.deepcopy(config), 'config': copy.deepcopy(config),
                'build': record(source), 'implementation': old_code, 'selection': [],
                'files': [record(path / n) for n in ('analysis.yaml', 'pipeline.yaml')]}
    write_json(path / 'analysis.json', manifest)
    write_json(path / 'checksums.json', {'analysis.json': record(path / 'analysis.json')['sha256']})
    monkeypatch.setattr(module, 'load', lambda p: json.loads((p / 'analysis.json').read_text()))
    monkeypatch.setattr(module, 'load_complete', lambda p: {})
    monkeypatch.setattr(module, 'verified_imports', lambda *args: [])
    def index(completion, destination, completed, species):
        write_json(destination / 'index.json', {'source': str(completion)})
    monkeypatch.setattr(module, 'create_index', index)
    original_check_output = subprocess.check_output
    monkeypatch.setattr(subprocess, 'check_output', lambda command, **kwargs:
        f'123|{root}\n' if command[0] == 'squeue' else original_check_output(command, **kwargs))
    monkeypatch.setenv('SLURM_JOB_ID', '123')  # The allocated audit itself is allowed.
    preview = module.upgrade(path, min_sample_fraction=requested)
    assert preview['state'] == 'ready'
    assert not list((path / 'jobs').glob('code_migration_*'))
    assert (path / 'pipeline.yaml').read_bytes() == previous_config
    result = module.upgrade(path, apply=True, min_sample_fraction=requested)
    assert result['state'] == 'upgraded'
    expected = 0 if requested is None else requested
    current = json.loads((path / 'analysis.json').read_text())
    for section, filename in [('pipeline', 'pipeline.yaml'), ('config', 'analysis.yaml')]:
        assert current[section]['alignment'] == {'enabled': True, 'min_sample_fraction': expected}
        assert current[section]['unchanged'] == config['unchanged']
        assert yaml.safe_load((path / filename).read_text()) == current[section]
    for entry in current['files']:
        verify(entry)
    assert json.loads((path / 'checksums.json').read_text())['analysis.json'] == record(path / 'analysis.json')['sha256']
    assert (path / 'jobs/code_migration_0003/pipeline.yaml').read_bytes() == previous_config
    assert current['code_migrations'][-1]['analysis_conditions'] == {'before': 0, 'after': expected}
    assert record(output) == old_output
    assert module.upgrade(path, min_sample_fraction=requested)['state'] == 'already_upgraded'


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
