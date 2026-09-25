"""Verify the reusable-build boundary with real Snakemake and synthetic tools."""
import copy
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

import analysis
import dataset
from build_products import complete, load_complete
from common import read_tsv, write_tsv
from test_datasets import dataset_project, imported, new_dataset


def test_analysis_requires_completed_build_and_rejects_build_settings(dataset_project):
    root = dataset_project
    imported(root)
    build = dataset.prepare(root,'incomplete',root/'config/build.yaml')
    with pytest.raises(ValueError,match='build is incomplete'):
        analysis.prepare(root,'test',root/'config/analysis.yaml',build)
    assert not (root/'analyses/test').exists()
    with pytest.raises(ValueError,match='reserved for build outputs'):
        analysis.prepare(root,'build_incomplete',root/'config/analysis.yaml',build)
    bad = root/'bad.yaml'; bad.write_text('odb:\n  node: 1\n')
    with pytest.raises(ValueError,match='unknown analysis settings'):
        analysis.settings(root,bad,build)
    with pytest.raises(ValueError,match='mapping incomplete'):
        complete(build)
    assert not (build/'completed.json').exists()


def test_retry_resources_do_not_change_frozen_build(dataset_project):
    root = dataset_project
    build = new_dataset(root)
    before = (build/'build.json').read_bytes()
    override = root/'retry.yaml'
    override.write_text(yaml.safe_dump({'translation':{'table':2},'slurm':{'stages':{'assembly':{'mem_mb':256000,'time':'7-00:00:00'}}}}))
    commands = dataset.submit(build,until='assembly',dry_run=True,resources=override)
    assert '--mem=256000M' in commands[0] and '--time=7-00:00:00' in commands[0]
    assert (build/'build.json').read_bytes() == before
    assert dataset.load(build)['analysis']['translation']['table'] == 1
    frozen = json.loads((build/'jobs/submission_0001.resources.json').read_text())
    assert frozen['slurm']['stages']['assembly']['mem_mb'] == 256000
    override.write_text('slurm:\n  stages:\n    assembly:\n      cpus: 0\n')
    with pytest.raises(ValueError,match='invalid assembly.cpus'):
        dataset.submit(build,until='assembly',dry_run=True,resources=override)


def test_new_native_products_reject_changed_biological_conditions(dataset_project):
    root = dataset_project
    build = new_dataset(root)
    dataset.submit(build,until='quant',dry_run=True)
    for step in dataset.STAGES: dataset.worker(build,step,1)
    config = root/'config/build.yaml'
    cfg = yaml.safe_load(config.read_text()); cfg['genegalleon']['settings']['assembly_method'] = 'Trinity'
    config.write_text(yaml.safe_dump(cfg))
    state = dataset.plan(root,config,'input/new.tsv')[-1][0]
    assert state['assembly'] == 'conflict'
    assert 'assembly settings differ' in state['reason']


def test_complete_build_multiple_analyses_and_species_updates(dataset_project,fake_odb,frozen_reference,command_environment,monkeypatch):
    root = dataset_project
    snakemake = shutil.which('snakemake'); seqkit = shutil.which('seqkit')
    if not snakemake or not seqkit: pytest.skip('Snakemake and seqkit required')
    imported(root)
    taxonomy = root/'resources/taxonomy/taxa.sqlite'; taxonomy.parent.mkdir(parents=True)
    shutil.copy2(root/'input/taxa.sqlite',taxonomy)
    reference = root/'resources/orthodb/v12_3193'; reference.parent.mkdir(parents=True)
    reference.symlink_to(frozen_reference,target_is_directory=True)
    events = root/'mapping_events.txt'
    env = command_environment({'python':sys.executable,'seqkit':seqkit,'ODB-mapper':fake_odb})
    env['FAKE_ODB_LOG'] = str(events)
    def execute(path,target, mapping=False):
        cmd = [snakemake,'--snakefile',str(root/'workflow/Snakefile'),'--configfile',str(path/'pipeline.yaml'),
               '--cores','2','--resources','mem_mb=16000']
        if mapping: cmd += ['--set-threads','odb_map=1','--set-resources','odb_map:mem_mb=3000']
        result = subprocess.run([*cmd,'--',target],cwd=root,env=env,capture_output=True,text=True,timeout=120)
        logs = '\n'.join(str(p)+': '+p.read_text()[-2000:] for p in (root/'logs').rglob('*.log'))
        assert result.returncode == 0, result.stdout+result.stderr+logs
    build = dataset.prepare(root,'base',root/'config/build.yaml')
    dataset.materialize(build)
    execute(build,'mapping',mapping=True)
    receipt = complete(build)
    original = receipt.read_bytes()
    assert dataset.submit(build,dry_run=True) == []
    assert set(load_complete(build)['products']) == {'Alpha_plant','Beta_sp-X','Gamma_plant'}
    assert len(events.read_text().splitlines()) == 1
    # Both analyses use the same completed proteins/mappings, with no mapper or translator rules.
    cfg = yaml.safe_load((root/'config/analysis.yaml').read_text())
    cfg['inputs']['species_trait'] = None
    cfg['phylogeny']['trees'] = []; cfg['phylogeny']['contrast_pairs']['enabled'] = False
    config = root/'analysis.local.yaml'
    for name,threshold,expected in [('loose',0.5,{'Alpha_plant','Beta_sp-X'}),('strict',0.7,{'Alpha_plant'})]:
        cfg['selection']['busco_threshold'] = threshold
        config.write_text(yaml.safe_dump(cfg))
        run = analysis.prepare(root,name,config,build)
        execute(run,'all')
        execute(run,'phenoradar_inputs')
        output = root/'results'/name
        assert {r['species'] for r in read_tsv(output/'orthogroups/expression/tpm.tsv')} == expected
        assert {r['species'] for r in read_tsv(output/'phenoradar_inputs/tpm.tsv')} == expected
        qc = json.loads((output/'orthogroups/mapping/merge_qc.json').read_text())
        assert qc['mode'] == 'completed_build' and qc['mapped_species'] == 0
        assert len(events.read_text().splitlines()) == 1
        assert receipt.read_bytes() == original
    # Resource changes are recorded per analysis submission without editing the frozen analysis.
    frozen = (run/'analysis.json').read_bytes()
    override = root/'resources.yaml'; override.write_text('slurm:\n  stages:\n    controller:\n      time: "5-00:00:00"\n  rules:\n    infer_species_tree:\n      mem_mb: 64000\n')
    cmd = analysis.submit(run,dry_run=True,resources=override)
    assert '--time=5-00:00:00' in cmd
    assert (run/'analysis.json').read_bytes() == frozen
    monkeypatch.setattr(subprocess,'check_output',lambda cmd,**kwargs: '123;cluster\n' if cmd[0]=='sbatch' else '123\n')
    analysis.submit(run,resources=override)
    with pytest.raises(ValueError,match='queued/running'):
        analysis.submit(run)
    # Removing then restoring species creates new build membership while preserving cached mappings.
    metadata = root/'input/metadata.tsv'; original_rows = read_tsv(metadata)
    for name,rows in [('removed',[original_rows[0]]),('restored',original_rows)]:
        policy = root/'config/excluded_accessions.tsv'
        blocked = [r for r in original_rows if r not in rows]
        write_tsv(policy,['accession','reason'],[{'accession':r['run'],'reason':'manual_failure_decision'} for r in blocked])
        updated = dataset.prepare(root,name,root/'config/build.yaml')
        assert dataset.submit(updated,until='quant',dry_run=True) == []
        dataset.materialize(updated); execute(updated,'mapping',mapping=True); complete(updated)
        assert set(load_complete(updated)['products']) == {r['scientific_name'].replace(' ','_') for r in rows}
        assert len(events.read_text().splitlines()) == 1
        assert [r['run'] for r in load_complete(updated)['excluded_runs']] == [r['run'] for r in blocked]
        if name == 'removed':
            subset = analysis.prepare(root,'without_excluded',config,updated)
            execute(subset,'all'); execute(subset,'phenoradar_inputs')
            assert {r['species'] for r in read_tsv(root/'results/without_excluded/phenoradar_inputs/tpm.tsv')} == {'Alpha_plant'}
            assert len(events.read_text().splitlines()) == 1
    # A completion record never legitimizes changed or deleted artifacts.
    protein = Path(load_complete(build)['products']['Alpha_plant']['protein']['path'])
    protein.write_text('changed\n')
    with pytest.raises(ValueError,match='registered file changed'):
        analysis.load(run)


def test_migration_preserves_caches_scientific_settings_and_never_overwrites(dataset_project):
    from migrate_phase_config import migrate
    root = dataset_project
    cfg = yaml.safe_load((root/'config/build.yaml').read_text())
    old = {k:copy.deepcopy(cfg[k]) for k in ('metadata','store','genegalleon')}
    old['odb_chunk_size'] = 7
    slurm = copy.deepcopy(cfg['slurm'])
    slurm['stages']['downstream'] = slurm['stages'].pop('controller')
    slurm['downstream_jobs'] = slurm.pop('jobs')
    slurm.pop('rules'); slurm.pop('default_resources')
    old['slurm'] = slurm
    old_path = root/'old-dataset.yaml'; old_path.write_text(yaml.safe_dump(old))
    flat = yaml.safe_load((root/'workflow/pipeline_defaults.yaml').read_text())
    flat['selection']['busco_threshold'] = 0.75
    flat['translation']['table'] = 4
    flat['odb']['existing_results'] = 'resources/odb_existing/tlight'
    old_flat = root/'old-config.yaml'; old_flat.write_text(yaml.safe_dump(flat))
    outputs = migrate(root,old_path,old_flat)
    build, lower = [yaml.safe_load(p.read_text()) for p in outputs]
    assert build['store'] == cfg['store']
    assert build['translation']['table'] == 4
    assert build['odb']['chunk_size'] == 7
    assert build['odb']['existing_results'] == 'resources/odb_existing/tlight'
    assert lower['selection']['busco_threshold'] == 0.75
    assert 'lineage' not in lower['phylogeny'] and 'translation' not in lower
    previous = [p.read_bytes() for p in outputs]
    with pytest.raises(ValueError,match='existing files are not overwritten'):
        migrate(root,old_path,old_flat)
    assert [p.read_bytes() for p in outputs] == previous


def test_submission_profile_applies_resources_to_rule_jobs(tmp_path):
    from phase_config import write_profile
    import yaml
    slurm = {'jobs':3, 'partition':'compute', 'account':'test', 'default_resources':{'mem_mb':4000,'runtime':60},
             'rules':{'odb_map':{'cpus':8,'mem_mb':64000,'runtime':120}}}
    profile = write_profile(tmp_path/'profile',slurm)
    cfg = yaml.safe_load((profile/'config.yaml').read_text())
    assert cfg['set-threads']['odb_map'] == 8
    assert cfg['set-resources']['odb_map'] == {'mem_mb':64000,'runtime':120}
    assert cfg['default-resources']['slurm_partition'] == 'compute'
    assert cfg['default-resources']['slurm_account'] == 'test'
