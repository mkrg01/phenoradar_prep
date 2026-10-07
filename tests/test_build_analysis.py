"""Verify the reusable-build boundary with real Snakemake and synthetic tools."""
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


@pytest.mark.parametrize('mode', ['species', 'sample', 'null', 'omitted'])
def test_analysis_species_list_is_enabled_by_input_path(dataset_project, mode):
    root = dataset_project
    build = imported(root)
    config = root / 'config/analysis.yaml'
    cfg = yaml.safe_load(config.read_text())
    assert 'species_list' not in cfg['selection']
    cfg['inputs']['species_trait'] = None
    subset = root / 'input/species_list.txt'
    subset.write_text('Beta_sp-X_B1\nGamma_plant_G1\n' if mode == 'sample'
                      else 'Beta_sp-X\nGamma_plant\n')
    enabled = mode in {'species', 'sample'}
    if mode == 'omitted':
        cfg['inputs'].pop('species_list')
    else:
        cfg['inputs']['species_list'] = 'input/species_list.txt' if enabled else None
    config.write_text(yaml.safe_dump(cfg))

    _, resolved, _, _, requested, report = analysis.settings(root, config, build)
    assert resolved['selection']['species_list'] is enabled
    expected = {'Beta_sp-X_B1', 'Gamma_plant_G1'}
    if not enabled:
        expected.add('Alpha_plant_A1')
    assert set(requested) == expected
    assert {r['species'] for r in report if r['selected']} == expected - {'Gamma_plant_G1'}
    reasons = {r['species']: r['reason'] for r in report}
    assert reasons['Gamma_plant_G1'] == 'below_busco_threshold'
    assert reasons['Alpha_plant_A1'] == ('outside_species_list' if enabled else '')

    run = analysis.prepare(root, mode, config, build)
    assert set((run / 'input/species_list.txt').read_text().splitlines()) == expected
    assert analysis.load(run)['selection'] == report


def test_analysis_requires_completed_build_and_rejects_build_settings(dataset_project):
    root = dataset_project
    imported(root)
    build = dataset.prepare(root,'incomplete',root/'config/build.yaml')
    with pytest.raises(ValueError,match='build is incomplete'):
        analysis.prepare(root,'test',root/'config/analysis.yaml',build)
    assert not (build/'downstream').exists()
    bad = root/'bad.yaml'; bad.write_text('odb:\n  ncbi_tax_id: 1\n')
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
    override.write_text(yaml.safe_dump({'translation':{'table':2},'slurm':{'stages':{'sample':{'mem_gb':256,'time':'7-00:00:00'}}}}))
    commands = dataset.submit(build,until='quant',dry_run=True,resources=override)
    assert '--mem=256000M' in commands[0] and '--time=7-00:00:00' in commands[0]
    assert (build/'build.json').read_bytes() == before
    assert dataset.load(build)['analysis']['translation']['table'] == 1
    frozen = json.loads((build/'jobs/submission_0001.resources.json').read_text())
    assert frozen['slurm']['stages']['sample']['mem_gb'] == 256
    for key, value in [('cpus', 0), ('mem_gb', 0), ('mem_gb', True), ('mem_gb', 1.5)]:
        override.write_text(yaml.safe_dump({'slurm': {'stages': {'sample': {key: value}}}}))
        with pytest.raises(ValueError,match=f'invalid sample.{key}'):
            dataset.submit(build,until='quant',dry_run=True,resources=override)


def test_new_native_products_reject_changed_biological_conditions(dataset_project):
    root = dataset_project
    build = new_dataset(root)
    dataset.submit(build,until='quant',dry_run=True)
    dataset.worker(build, 1)
    config = root/'config/build.yaml'
    from database_fixtures import database_from_stages
    frozen = dataset.load(build)
    source = database_from_stages(root, frozen['config']['store'], frozen['items'], 'conditions_source')
    cfg = yaml.safe_load(config.read_text()); cfg['reuse_from'] = str(source)
    cfg['genegalleon']['settings']['assembly_method'] = 'Trinity'
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
        logs = '\n'.join(str(p)+': '+p.read_text()[-2000:] for p in (root/'results').rglob('*.log'))
        assert result.returncode == 0, result.stdout+result.stderr+logs
    build = dataset.prepare(root,'base',root/'config/build.yaml')
    dataset.materialize(build)
    execute(build,'mapping',mapping=True)
    with pytest.raises(ValueError, match='expression incomplete'):
        complete(build)
    execute(build,'database',mapping=True)
    assert (build/'work/database/orthogroups/expression/runs/A1.tsv').is_file()
    assert not (root/'builds').exists()
    receipt = complete(build)
    original = receipt.read_bytes()
    assert dataset.submit(build,dry_run=True) == []
    assert set(load_complete(build)['products']) == {'Alpha_plant_A1','Beta_sp-X_B1','Gamma_plant_G1'}
    assert not events.exists()
    # Both analyses use the same completed proteins/mappings, with no mapper or translator rules.
    cfg = yaml.safe_load((root/'config/analysis.yaml').read_text())
    cfg['inputs']['species_trait'] = None
    cfg['phylogeny']['trees'] = []; cfg['phylogeny']['contrast_pairs']['enabled'] = False
    config = root/'analysis.local.yaml'
    for name,threshold,expected in [('loose',0.5,{'Alpha_plant_A1','Beta_sp-X_B1'}),('strict',0.7,{'Alpha_plant_A1'})]:
        cfg['selection']['busco_threshold'] = threshold
        config.write_text(yaml.safe_dump(cfg))
        run = analysis.prepare(root,name,config,build)
        execute(run,'all')
        execute(run,'phenoradar_inputs')
        assert run == build/'downstream'/name
        output = run
        assert 'tpm' not in analysis.load(run)['pipeline']
        assert load_complete(build)['tpm'] == {'multimap': 'error'}
        assert (output/'orthogroups/expression/runs/A1.tsv').samefile(build/'database/expression/runs/A1.tsv')
        assert {r['species'] for r in read_tsv(output/'orthogroups/expression/tpm.tsv')} == expected
        assert {r['species'] for r in read_tsv(output/'phenoradar_inputs/tpm.tsv')} == expected
        # Expression-only downstream runs do not even stage proteins or mappings.
        assert not (output/'orthogroups/mapping').exists()
        assert not (output/'proteins').exists()
        execute(run,'mapping')
        qc = json.loads((output/'orthogroups/mapping/snapshot.json').read_text())['qc']
        assert qc['mode'] == 'existing' and qc['mapped_species'] == 0
        assert not events.exists()
        assert receipt.read_bytes() == original
    # Resource changes are recorded per analysis submission without editing the frozen analysis.
    frozen = (run/'analysis.json').read_bytes()
    override = root/'resources.yaml'; override.write_text('slurm:\n  stages:\n    controller:\n      time: "5-00:00:00"\n  rules:\n    infer_species_tree:\n      mem_gb: 64\n')
    cmd = analysis.submit(run,dry_run=True,resources=override)
    assert '--time=5-00:00:00' in cmd
    assert (run/'analysis.json').read_bytes() == frozen
    monkeypatch.setattr(subprocess,'check_output',lambda cmd,**kwargs: '123;cluster\n' if cmd[0]=='sbatch' else '123\n')
    analysis.submit(run,resources=override)
    with pytest.raises(ValueError,match='queued/running'):
        analysis.submit(run)
    build_config = root/'config/build.yaml'
    value = yaml.safe_load(build_config.read_text()); value['reuse_from'] = str(build/'database')
    build_config.write_text(yaml.safe_dump(value))
    # Removing then restoring species creates new build membership while preserving cached mappings.
    metadata = root/'input/metadata.tsv'; original_rows = read_tsv(metadata)
    for name,rows in [('removed',[original_rows[0]]),('restored',original_rows)]:
        policy = root/'config/excluded_accessions.tsv'
        blocked = [r for r in original_rows if r not in rows]
        write_tsv(policy,['accession','reason'],[{'accession':r['run'],'reason':'manual_failure_decision'} for r in blocked])
        updated = dataset.prepare(root,name,root/'config/build.yaml')
        assert dataset.submit(updated,until='quant',dry_run=True) == []
        dataset.materialize(updated); execute(updated,'database',mapping=True); complete(updated)
        assert set(load_complete(updated)['products']) == {r['scientific_name'].replace(' ','_')+'_'+r['run'] for r in rows}
        assert not events.exists()
        assert [r['run'] for r in load_complete(updated)['excluded_runs']] == [r['run'] for r in blocked]
        if name == 'removed':
            subset = analysis.prepare(root,'without_excluded',config,updated)
            execute(subset,'all'); execute(subset,'phenoradar_inputs')
            assert {r['species'] for r in read_tsv(subset/'phenoradar_inputs/tpm.tsv')} == {'Alpha_plant_A1'}
            assert not events.exists()
        assert (updated/'database/expression/runs/A1.tsv').samefile(build/'database/expression/runs/A1.tsv')
        assert analysis.prepare(root,'loose',config,updated) == updated/'downstream/loose'
    assert not (root/'analyses').exists()
    assert not (root/'work').exists()
    assert not (root/'logs').exists()
    # A completion record never legitimizes changed or deleted artifacts.
    protein = Path(load_complete(build)['products']['Alpha_plant_A1']['protein']['path'])
    protein.write_text('changed\n')
    with pytest.raises(ValueError,match='registered file changed'):
        analysis.load(run)


def test_submission_profile_applies_resources_to_rule_jobs(tmp_path):
    from phase_config import write_profile
    import yaml
    slurm = {'jobs':3, 'partition':'compute', 'default_resources':{'mem_gb':4,'time':'01:00:00'},
             'rules':{'odb_map':{'cpus':8,'mem_gb':64,'time':'02:00:00'}}}
    profile = write_profile(tmp_path/'profile',slurm)
    cfg = yaml.safe_load((profile/'config.yaml').read_text())
    assert cfg['set-threads']['odb_map'] == 8
    assert cfg['set-resources']['odb_map'] == {'mem_mb':64000,'runtime':120}
    assert cfg['default-resources'] == {'mem_mb':4000, 'runtime':60, 'slurm_partition':'compute'}
    assert cfg['slurm-no-account'] is True
    # Conversion must not alter the saved input units or compound on reuse.
    assert slurm['default_resources']['mem_gb'] == 4
    assert slurm['rules']['odb_map']['mem_gb'] == 64
    assert slurm['default_resources']['time'] == '01:00:00'
    assert slurm['rules']['odb_map']['time'] == '02:00:00'
    write_profile(profile, slurm)
    assert yaml.safe_load((profile/'config.yaml').read_text()) == cfg
