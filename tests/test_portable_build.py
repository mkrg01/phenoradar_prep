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
from incremental_odb import import_snapshot, plan
from portable_build import completion_path, register_products
from test_datasets import dataset_project, imported
from test_incremental_odb import snapshot


def test_register_odb_cache_is_automatic_idempotent_and_independent(tmp_path):
    proteins = tmp_path / 'proteins'; proteins.mkdir()
    for species in ('Alpha_plant', 'Beta_plant'):
        (proteins / f'{species}_protein.fa').write_text(f'>{species}_g1\nMK\n>{species}_g2\nMP\n')
    source = tmp_path / 'old'
    snapshot(source, proteins, ['Alpha_plant'])
    cache = tmp_path / 'cache'
    registered = import_snapshot(source, cache)
    assert import_snapshot(source, cache) == registered
    # Historical path/timestamps are not mapping identities.
    d = json.loads((source / 'snapshot.json').read_text())
    d['proteins'][0]['path'] = '/retired/protein.fa'
    d['created_at'] = 'later import of the same result'
    write_json(source / 'snapshot.json', d)
    assert import_snapshot(source, cache) == registered
    shutil.rmtree(source)
    samples = tmp_path / 'samples.tsv'
    write_tsv(samples, ['species','odb_species'], [{'species':s,'odb_species':s} for s in ['Alpha_plant','Beta_plant']])
    result = plan(samples, proteins, tmp_path / 'plan', cache)
    assert result['reused_species'] == ['Alpha_plant']
    assert result['mapped_species'] == ['Beta_plant']
    assert len(list(cache.glob('*/*/snapshot.json'))) == 1
    (registered / 'annotations.tsv').write_text('corrupted\n')
    with pytest.raises(ValueError, match='changed after completion'):
        plan(samples, proteins, tmp_path / 'bad-plan', cache)


def test_register_rejects_wrong_odb_node_and_corrupt_source(tmp_path):
    proteins = tmp_path / 'proteins'; proteins.mkdir()
    (proteins / 'Alpha_plant_protein.fa').write_text('>Alpha_plant_g1\nMK\n')
    source = tmp_path / 'source'; snapshot(source, proteins, ['Alpha_plant'])
    with pytest.raises(ValueError, match='version/node differs'):
        import_snapshot(source, tmp_path / 'cache', node=33090)
    (source / 'annotations.tsv').write_text('changed\n')
    with pytest.raises(ValueError, match='changed after completion'):
        import_snapshot(source, tmp_path / 'cache')
    assert not list((tmp_path / 'cache').glob('*/*/snapshot.json'))


def execute(root, path, target, env, mapping=False):
    cmd = [shutil.which('snakemake'), '--snakefile', str(root/'workflow/Snakefile'),
           '--configfile', str(path/'pipeline.yaml'), '--cores','2','--resources','mem_mb=16000']
    if mapping: cmd += ['--set-threads','odb_map=1','--set-resources','odb_map:mem_mb=3000']
    result = subprocess.run([*cmd,'--',target],cwd=root,env=env,capture_output=True,text=True,timeout=120)
    if result.returncode:
        logs = '\n'.join(str(p)+': '+p.read_text()[-1500:] for p in (root/'logs').rglob('*.log'))
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
    execute(root,build,'mapping',env,mapping=True)
    complete(build)
    return root,build,env,events


def test_bundle_relocation_analysis_and_incremental_reuse(completed_project):
    root,build,env,events = completed_project
    source = build/'products'
    assert (source/'manifest.json').exists()
    assert complete(build) == build/'completed.json'
    # A crash after products publication but before the pointer is recoverable.
    (build/'completed.json').unlink()
    assert complete(build) == build/'completed.json'
    original = json.loads((source/'manifest.json').read_text())
    assert original['schema_version'] == 3
    odb_snapshot = json.loads((source/'odb/snapshot.json').read_text())
    assert len(odb_snapshot['reference_sha256s']) == 1
    assert all(not Path(r['path']).is_absolute() for r in original['files'])
    second = root.parent/'second'; second.mkdir()
    for directory in ['workflow','config']:
        shutil.copytree(root/directory,second/directory)
    shutil.copy2(root/'run_pipeline.sh',second/'run_pipeline.sh')
    # Only products (and an independent taxonomy reference) are transported.
    moved = second/'builds/copied/products'
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
        execute(second,run,'all',env)
        execute(second,run,'phenoradar_inputs',env)
        assert {r['species'] for r in read_tsv(second/'results/copied/phenoradar_inputs/tpm.tsv')} == {'Alpha_plant','Beta_sp-X'}
        assert len(events.read_text().splitlines()) == 1
        cfg = yaml.safe_load((second/'config/build.yaml').read_text())
        cfg['metadata'] = str(moved/'metadata.tsv')
        (second/'config/build.yaml').write_text(yaml.safe_dump(cfg))
        rows,odb = register_products(moved,second/cfg['store'],second/cfg['odb']['cache_dir'],
                                    cfg['busco']['lineage'],cfg['translation'],cfg['odb']['node'])
        assert all(r['status']=='registered' for r in rows)
        assert odb.is_relative_to(second/'resources/odb_cache')
        changed_reference = second/'changed-reference.json'; changed_reference.write_text('changed')
        # Known native reference constraints survive transport and re-registration.
        copied_samples = second/'copied-samples.tsv'
        write_tsv(copied_samples, ['species','odb_species'],
                  [{'species':s,'odb_species':p['odb_species']} for s,p in loaded['products'].items()])
        with pytest.raises(ValueError,match='reference changed'):
            plan(copied_samples,moved/'proteins',second/'rejected-plan',second/cfg['odb']['cache_dir'],
                 reference=changed_reference)
        report = dataset.plan(second,second/'config/build.yaml')[-1]
        assert all(r[s]=='reuse' for r in report for s in dataset.STAGES)
        expansion = dataset.prepare(second,'expanded',second/'config/build.yaml')
        assert dataset.submit(expansion,until='quant',dry_run=True) == []
        dataset.materialize(expansion)
        execute(second,expansion,'mapping',env,mapping=True)
        translated = second/'results'/('build_' + expansion.name)/'proteins'
        assert all(json.loads(p.read_text())['reused'] for p in translated.glob('*_protein.json'))
        complete(expansion)
        assert len(events.read_text().splitlines()) == 1
        assert (expansion/'products/manifest.json').exists()
    finally:
        offline.rename(root)


def test_bundle_rejects_missing_corrupt_and_escaping_paths(completed_project):
    root,build,_,_ = completed_project
    dest = root.parent/'copied-products'; shutil.copytree(build/'products',dest)
    manifest = dest/'manifest.json'; original = json.loads(manifest.read_text())
    entry = original['products']['Alpha_plant']['protein']
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


def test_register_odb_only_cli_does_not_require_species_inputs(dataset_project):
    root = dataset_project
    proteins = root/'proteins'; proteins.mkdir()
    (proteins/'Alpha_plant_protein.fa').write_text('>Alpha_plant_g1\nMK\n')
    source = root/'old'; snapshot(source,proteins,['Alpha_plant'])
    shutil.rmtree(root/'input')
    command = [os.sys.executable,str(root/'workflow/scripts/dataset.py'),'register','--root',str(root)]
    result = subprocess.run([*command,'--odb-only','--odb-results','old'],capture_output=True,text=True)
    assert result.returncode == 0,result.stderr
    output = json.loads(result.stdout)
    assert output['registered'] == 0 and len(output['odb_snapshots']) == 1
    assert not (root/'resources/dataset_assets').exists()
    assert (Path(output['odb_snapshots'][0])/'snapshot.json').exists()
    invalid = subprocess.run([*command,'--odb-only'],capture_output=True,text=True)
    assert invalid.returncode != 0 and '--odb-only requires --odb-results' in invalid.stderr


def test_portable_conditions_and_settings_are_preserved(completed_project):
    root,build,_,_ = completed_project
    bundle = root.parent/'conditions'; shutil.copytree(build/'products',bundle)
    manifest = bundle/'manifest.json'; data = json.loads(manifest.read_text())
    cfg, _ = dataset.settings(root,root/'config/build.yaml')
    data['products']['Alpha_plant']['conditions'] = cfg['conditions']
    data['sha256'] = digest({k:v for k,v in data.items() if k != 'sha256'})
    write_json(manifest,data)
    store = root/'resources/copied-store'; cache = root/'resources/copied-odb'
    with pytest.raises(ValueError,match='settings differ'):
        register_products(bundle,store,cache,'other_odb12',cfg['translation'],cfg['odb']['node'])
    assert not store.exists()
    register_products(bundle,store,cache,cfg['busco']['lineage'],cfg['translation'],cfg['odb']['node'])
    config = yaml.safe_load((root/'config/build.yaml').read_text());config['store'] = str(store)
    config['metadata'] = str(bundle/'metadata.tsv')
    new_config = root/'copied.yaml';new_config.write_text(yaml.safe_dump(config))
    assert all(r['assembly']=='reuse' for r in dataset.plan(root,new_config)[-1])
    config['genegalleon']['settings']['assembly_method'] = 'changed-assembler'
    new_config.write_text(yaml.safe_dump(config))
    alpha = next(r for r in dataset.plan(root,new_config)[-1] if r['species']=='Alpha_plant')
    assert alpha['assembly']=='conflict' and 'assembly settings differ' in alpha['reason']
