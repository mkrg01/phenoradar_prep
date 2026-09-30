"""Copied completed builds work without access to their original project."""
import copy
import json
import os
import shutil
import subprocess
from pathlib import Path

import pytest
import yaml

import analysis
import dataset
from build_products import complete, load_complete
from common import read_tsv, write_json, write_tsv
from dataset_assets import digest
from portable_build import completion_path
from test_datasets import dataset_project, imported


def execute(root, path, target, env, mapping=False):
    cmd = [shutil.which('snakemake'), '--snakefile', str(root/'workflow/Snakefile'),
           '--configfile', str(path/'pipeline.yaml'), '--cores','2','--resources','mem_mb=16000']
    if mapping: cmd += ['--set-threads','odb_map=1','--set-resources','odb_map:mem_mb=3000']
    result = subprocess.run([*cmd,'--',target],cwd=root,env=env,capture_output=True,text=True,timeout=120)
    if result.returncode:
        logs = '\n'.join(str(p)+': '+p.read_text()[-1500:] for p in (root/'results').rglob('*.log'))
        pytest.fail(result.stdout+result.stderr+logs)


@pytest.fixture
def completed_project(dataset_project, fake_odb, frozen_reference, command_environment):
    root = dataset_project
    if not shutil.which('snakemake') or not shutil.which('seqkit'):
        pytest.skip('Snakemake and seqkit required')
    imported(root)
    taxonomy = root/'resources/taxonomy/taxa.sqlite'; taxonomy.parent.mkdir(parents=True)
    shutil.copy2(root/'input/taxa.sqlite', taxonomy)
    reference = root/'resources/orthodb/v12_3193'; reference.parent.mkdir(parents=True)
    reference.symlink_to(frozen_reference, target_is_directory=True)
    env = command_environment({'python':os.sys.executable,'seqkit':shutil.which('seqkit'),'ODB-mapper':fake_odb})
    events = root.parent/'events.txt'; env['FAKE_ODB_LOG'] = str(events)
    build = dataset.prepare(root,'base',root/'config/build.yaml')
    dataset.materialize(build)
    execute(root,build,'database',env,mapping=True)
    complete(build)
    return root,build,env,events


def test_bundle_relocation_and_analysis(completed_project):
    root,build,env,events = completed_project
    source = build/'database'
    assert (source/'manifest.json').exists()
    assert complete(build) == build/'completed.json'
    # A crash after products publication but before the pointer is recoverable.
    (build/'completed.json').unlink()
    assert complete(build) == build/'completed.json'
    original = json.loads((source/'manifest.json').read_text())
    assert original['schema_version'] == 5
    assert original['tpm'] == {'multimap': 'error'}
    assert (source/'expression/runs/A1.tsv').is_file()
    odb_snapshot = json.loads((source/'odb/snapshot.json').read_text())
    assert odb_snapshot['reference_sha256s'] == []
    assert all(not Path(r['path']).is_absolute() for r in original['files'])
    second = root.parent/'second'; second.mkdir()
    for directory in ['workflow','config']:
        shutil.copytree(root/directory,second/directory)
    shutil.copy2(root/'run_pipeline.sh',second/'run_pipeline.sh')
    # Only products (and an independent taxonomy reference) are transported.
    moved = second/'results/copied/database'
    shutil.copytree(source,moved)
    taxonomy = second/'resources/taxonomy'; taxonomy.mkdir(parents=True)
    shutil.copy2(root/'resources/taxonomy/taxa.sqlite',taxonomy/'taxa.sqlite')
    cfg = yaml.safe_load((second/'config/analysis.yaml').read_text())
    cfg['inputs']['species_trait'] = None
    cfg['phylogeny']['trees'] = []; cfg['phylogeny']['contrast_pairs']['enabled'] = False
    (second/'config/analysis.yaml').write_text(yaml.safe_dump(cfg))
    offline = root.with_name('original-unavailable')
    root.rename(offline)
    try:
        loaded = load_complete(moved)
        assert all(Path(p['path']).is_relative_to(moved) for p in loaded['files'])
        assert completion_path(moved.parent) == moved/'manifest.json'
        run = analysis.prepare(second,'copied',second/'config/analysis.yaml',moved)
        assert run == second/'results/copied/downstream/copied'
        execute(second,run,'all',env)
        execute(second,run,'phenoradar_inputs',env)
        assert {r['species'] for r in read_tsv(run/'phenoradar_inputs/tpm.tsv')} == {'Alpha_plant_A1','Beta_sp-X_B1'}
        assert not (run/'orthogroups/expression/runs/A1.qc.json').samefile(moved/'expression/runs/A1.qc.json')
        assert not events.exists()
    finally:
        offline.rename(root)


def test_bundle_rejects_missing_corrupt_and_escaping_paths(completed_project):
    root,build,_,_ = completed_project
    dest = root.parent/'copied-products'; shutil.copytree(build/'database',dest)
    manifest = dest/'manifest.json'; original = json.loads(manifest.read_text())
    entry = original['products']['Alpha_plant_A1']['protein']
    protein = dest/entry['path']; original_bytes = protein.read_bytes()
    protein.unlink()
    with pytest.raises(ValueError,match='registered file missing'): load_complete(dest)
    protein.write_bytes(original_bytes+b'changed')
    with pytest.raises(ValueError,match='registered file changed'): load_complete(dest)
    protein.write_bytes(original_bytes)
    for replacement in ['../outside.fa','/tmp/outside.fa']:
        bad = copy.deepcopy(original)
        bad['files'][0]['path'] = replacement
        bad['sha256'] = digest({k:v for k,v in bad.items() if k!='sha256'})
        write_json(manifest,bad)
        with pytest.raises(ValueError,match='invalid relative product path'): load_complete(dest,verify_files=False)
    write_json(manifest,original)
    outside = root.parent/'outside.fa'; outside.write_bytes(original_bytes)
    protein.unlink(); protein.symlink_to(outside)
    with pytest.raises(ValueError,match='escapes bundle'): load_complete(dest,verify_files=False)


def test_legacy_products_remain_usable(completed_project):
    root, build, env, _ = completed_project
    legacy = root / 'builds/legacy/products'
    shutil.copytree(build / 'database', legacy)
    manifest = legacy / 'manifest.json'
    old = json.loads(manifest.read_text())
    old.update(schema_version=4, build_id='legacy')
    old.pop('tpm')
    old['files'] = [entry for entry in old['files'] if not entry['path'].startswith('expression/')]
    for product in old['products'].values():
        for key in ('expression', 'expression_qc', 'source_expression_qc_sha256'):
            product.pop(key)
    old['sha256'] = digest({k:v for k,v in old.items() if k != 'sha256'})
    write_json(manifest, old)
    shutil.rmtree(legacy / 'expression')
    config = root / 'legacy-analysis.yaml'
    config.write_text('inputs:\n  species_trait: null\nphylogeny:\n  trees: []\n  contrast_pairs:\n    enabled: false\n')
    assert load_complete(legacy.parent)['schema_version'] == 4
    downstream = analysis.prepare(root, 'legacy', config, legacy.parent)
    assert downstream == root / 'results/legacy/downstream/legacy'
    execute(root, downstream, 'all', env)
    assert (downstream / 'orthogroups/expression/tpm.tsv').is_file()


def test_database_expression_and_policy_are_required(completed_project):
    root, build, _, _ = completed_project
    config = root / 'bad-analysis.yaml'
    config.write_text('tpm:\n  multimap: split\n')
    with pytest.raises(ValueError, match='unknown analysis settings: tpm'):
        analysis.prepare(root, 'bad', config, build)
    database = root.parent / 'incomplete-database'
    shutil.copytree(build / 'database', database)
    manifest = database / 'manifest.json'
    original = json.loads(manifest.read_text())
    for policy in ('drop', 'split'):
        bad = copy.deepcopy(original)
        bad['tpm']['multimap'] = policy
        bad['sha256'] = digest({k:v for k,v in bad.items() if k != 'sha256'})
        write_json(manifest, bad)
        with pytest.raises(ValueError, match='database TPM policy must be error'):
            load_complete(database, verify_files=False)
    write_json(manifest, original)
    product = original['products']['Alpha_plant_A1']
    expression = database / product['expression']['path']
    expression.write_text('corrupt\n')
    with pytest.raises(ValueError, match='registered file changed'):
        load_complete(database)
    bad = copy.deepcopy(original)
    bad['products']['Alpha_plant_A1'].pop('expression')
    bad['sha256'] = digest({k:v for k,v in bad.items() if k != 'sha256'})
    write_json(manifest, bad)
    with pytest.raises(ValueError, match='incomplete database product'):
        load_complete(database, verify_files=False)


def test_downstream_completion_records_outputs_without_its_own_receipt(completed_project, monkeypatch):
    root, build, _, _ = completed_project
    config = root / 'completion-analysis.yaml'
    config.write_text('inputs:\n  species_trait: null\nphylogeny:\n  trees: []\n  contrast_pairs:\n    enabled: false\n')
    downstream = analysis.prepare(root, 'finished', config, build)
    def run_pipeline(command, **kwargs):
        if command[-1] == 'all':
            write_json(downstream / 'run.json', {'completed': True})
        else:
            assert command[-1] == 'phenoradar_inputs'
            write_tsv(downstream / 'phenoradar_inputs/tpm.tsv', ['species', 'orthogroup', 'tpm'],
                      [{'species':'Alpha_plant_A1', 'orthogroup':'OG1', 'tpm':1000000}])
    monkeypatch.setattr(analysis.subprocess, 'run', run_pipeline)
    analysis.run(downstream, local=True)
    assert analysis.status(downstream)['state'] == 'complete'
    receipt = json.loads((downstream / 'completed.json').read_text())
    paths = {Path(entry['path']).relative_to(downstream).as_posix() for entry in receipt['files']}
    assert paths == {'run.json', 'phenoradar_inputs/tpm.tsv'}
    analysis.run(downstream, local=True)
    assert analysis.status(downstream)['state'] == 'complete'
    (downstream / 'phenoradar_inputs/tpm.tsv').write_text('changed\n')
    with pytest.raises(ValueError, match='registered file changed'):
        analysis.status(downstream)
