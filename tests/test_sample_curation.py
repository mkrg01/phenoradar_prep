"""QC recording preserves successful work across exclusions and partial builds."""
import json
import os
from pathlib import Path
import shutil
import subprocess

import pytest
import yaml

import dataset
import metadata_catalog as catalog
from common import read_tsv, write_tsv
from test_datasets import dataset_project, fake_genegalleon, native_events


@pytest.fixture
def sample_project(dataset_project):
    root = dataset_project
    fake_genegalleon(root)
    rows = read_tsv(root / 'input/metadata.tsv')
    for row in rows:
        row.update(bioproject='P1', sample_group='leaf', exclusion='no', total_bases='100', lib_layout='paired')
    write_tsv(root / 'input/metadata.tsv', list(rows[0]), rows)
    folder = root / 'datasets/leaf'
    folder.mkdir(parents=True)
    (folder / 'select_rules.tsv').write_text('rule_id\tenabled\nleaf\tyes\n')
    definition = {'schema_version': 1, 'name': 'leaf', 'build_config': 'config/build.yaml',
                  'search_string': 'plants', 'sample_group': 'leaf', 'rule_set': 'plantae',
                  'rules': 'datasets/leaf/select_rules.tsv',
                  'previous_metadata': 'input/metadata.tsv', 'accepted_samples': 'datasets/leaf/accepted_samples.tsv',
                  'excluded_accessions': 'config/excluded_accessions.tsv'}
    config = folder / 'selection.yaml'
    config.write_text(yaml.safe_dump(definition))
    write_tsv(root / 'config/excluded_accessions.tsv', ['accession', 'reason'], [])
    cfg = catalog.configuration(root, config)
    build = dataset.prepare(root, 'first', root / 'config/build.yaml')
    dataset.submit(build, until='quant', dry_run=True)
    return root, cfg, build


def finish(build, index):
    dataset.worker(build, index)


def test_busco_failures_are_excluded_while_completed_passing_samples_are_adopted(sample_project, monkeypatch):
    root, cfg, build = sample_project
    monkeypatch.setenv('FAKE_GG_LOW_BUSCO', '1')
    dataset.worker(build, 1)
    monkeypatch.delenv('FAKE_GG_LOW_BUSCO')
    dataset.worker(build, 1)
    finish(build, 2)
    monkeypatch.setenv('FAKE_GG_FAIL_ASSEMBLY', '1')
    with pytest.raises(subprocess.CalledProcessError): dataset.worker(build, 3)
    before = (root / 'input/metadata.tsv').read_bytes()
    exclusions = (root / 'config/excluded_accessions.tsv').read_bytes()
    preview = catalog.record_successes(root, cfg, build, dry_run=True)
    assert preview['recorded_runs'] == ['B1']
    assert preview['excluded_runs'] == ['A1']
    assert not (root / cfg['accepted_samples']).exists()
    assert (root / cfg['excluded_accessions']).read_bytes() == exclusions
    assert preview['samples'][0]['busco_complete'] == 0
    assert preview['samples'][0]['status'] == 'busco_below_threshold'
    assert preview['samples'][2]['status'] == 'incomplete'
    assert preview['busco_threshold'] == 0.5
    assert not (build / 'completed.json').exists()
    result = catalog.record_successes(root, cfg, build)
    assert result['recorded_runs'] == ['B1']
    assert result['excluded_runs'] == ['A1']
    assert read_tsv(root / cfg['excluded_accessions']) == [
        dict(accession='A1', taxid=read_tsv(root / 'input/metadata.tsv')[0]['taxid'], bioproject='P1',
             reason='busco_completeness_below_0.5', source_build='first')]
    provenance = json.loads((root / cfg['accepted_samples']).with_suffix('.provenance.json').read_text())
    assert provenance['excluded_runs'] == ['A1']
    exclusion_bytes = (root / cfg['excluded_accessions']).read_bytes()
    repeated = catalog.record_successes(root, cfg, build, runs=['A1'])
    assert repeated['recorded_runs'] == repeated['excluded_runs'] == []
    assert repeated['samples'][0]['status'] == 'excluded'
    assert (root / cfg['excluded_accessions']).read_bytes() == exclusion_bytes
    assert [r['run'] for r in read_tsv(root / cfg['accepted_samples'])] == ['B1']
    assert (root / 'input/metadata.tsv').read_bytes() == before
    rows = read_tsv(root / 'input/metadata.tsv')
    larger = dict(rows[1], run='B_BIG', total_bases='9000')
    source = root / 'curated.tsv'
    write_tsv(source, list(rows[0]), [*rows, larger, dict(rows[0], run='A2', total_bases='50')])
    catalog.update_metadata(root, cfg, 'work/next_metadata', source)
    assert [r['run'] for r in read_tsv(root / 'input/metadata.tsv')] == ['A2', 'B1', 'G1']


def test_exact_threshold_passes_and_later_threshold_change_does_not_revoke_adoption(sample_project, monkeypatch):
    root, cfg, build = sample_project
    monkeypatch.setenv('FAKE_GG_HALF_BUSCO', '1')
    dataset.worker(build, 1)
    monkeypatch.delenv('FAKE_GG_HALF_BUSCO')
    dataset.worker(build, 1)
    result = catalog.record_successes(root, cfg, build, runs=['A1'])
    assert result['recorded_runs'] == ['A1']
    assert result['excluded_runs'] == []
    assert result['samples'][0]['busco_complete'] == 1
    assert result['samples'][0]['busco_total'] == 2
    evidence = (root / cfg['accepted_samples']).read_bytes()
    definition = yaml.safe_load(cfg['_path'].read_text())
    definition['busco_threshold'] = 0.75
    cfg['_path'].write_text(yaml.safe_dump(definition))
    cfg = catalog.configuration(root, cfg['_path'])
    result = catalog.record_successes(root, cfg, build, runs=['A1'])
    assert result['recorded_runs'] == []
    assert result['samples'][0]['status'] == 'busco_below_threshold'
    assert result['samples'][0]['retained_accepted'] is True
    assert result['excluded_runs'] == []
    assert read_tsv(root / cfg['excluded_accessions']) == []
    assert (root / cfg['accepted_samples']).read_bytes() == evidence
    rows = read_tsv(root / 'input/metadata.tsv')
    source = root / 'curated.tsv'
    write_tsv(source, list(rows[0]), [*rows, dict(rows[0], run='A_BIG', total_bases='9000')])
    catalog.update_metadata(root, cfg, 'work/raised_threshold', source)
    assert read_tsv(root / 'input/metadata.tsv')[0]['run'] == 'A1'
    write_tsv(root / cfg['excluded_accessions'], ['accession', 'reason'],
              [{'accession': 'A1', 'reason': 'misidentified'}])
    result = catalog.record_successes(root, cfg, build, runs=['A1'])
    assert result['samples'][0]['status'] == 'excluded'
    assert read_tsv(root / cfg['accepted_samples']) == []


@pytest.mark.parametrize('exclusion_mode', ['manual', 'busco'])
def test_exclusion_replaces_same_taxid_and_new_build_automatically_reuses_partial_successes(
        sample_project, monkeypatch, exclusion_mode):
    root, cfg, old = sample_project
    finish(old, 1)
    monkeypatch.setenv('FAKE_GG_FAIL_BUSCO', '1')
    with pytest.raises(subprocess.CalledProcessError): dataset.worker(old, 2)
    monkeypatch.delenv('FAKE_GG_FAIL_BUSCO')
    frozen = (old / 'build.json').read_bytes()
    if exclusion_mode == 'manual':
        write_tsv(root / 'config/excluded_accessions.tsv', ['accession', 'reason'],
                  [{'accession': 'G1', 'reason': 'unusable_reads'}])
    else:
        monkeypatch.setenv('FAKE_GG_LOW_BUSCO', '1')
        dataset.worker(old, 3)
        monkeypatch.delenv('FAKE_GG_LOW_BUSCO')
    result = catalog.record_successes(root, cfg, old)
    assert result['excluded_runs'] == (['G1'] if exclusion_mode == 'busco' else [])
    rows = read_tsv(root / 'input/metadata.tsv')
    replacement = dict(rows[2], run='G2', total_bases='50')
    source = root / 'curated.tsv'
    write_tsv(source, list(rows[0]), [*rows, replacement])
    catalog.update_metadata(root, cfg, 'work/replacement', source)
    assert [r['run'] for r in read_tsv(root / 'input/metadata.tsv')] == ['A1', 'B1', 'G2']
    config = root / 'config/build.yaml'
    definition = yaml.safe_load(config.read_text()); definition['reuse_from'] = 'auto'
    config.write_text(yaml.safe_dump(definition))
    report = dataset.plan(root, config, name='second')[-1]
    assert [r['assembly'] for r in report] == ['reuse', 'reuse', 'pending']
    assert [r['busco'] for r in report] == ['reuse', 'pending', 'pending']
    second = dataset.prepare(root, 'second', config)
    assert (old / 'build.json').read_bytes() == frozen
    old.rename(old.with_name('source_offline'))
    commands = dataset.submit(second, until='quant', dry_run=True)
    assert len(commands) == 1 and '--array=2,3' in commands[0]
    finish(second, 2)
    finish(second, 3)
    assert len(native_events(second)) == 2
    assert {r['run'] for r in read_tsv(dataset.materialize(second) / 'metadata.tsv')} == {'A1', 'B1', 'G2'}


@pytest.mark.parametrize('fields', [catalog.EXCLUSION_FIELDS + ['reviewer'], ['accession']])
def test_busco_exclusion_preserves_manual_decisions_and_custom_columns(sample_project, monkeypatch, fields):
    root, cfg, build = sample_project
    manual = {key: value for key, value in dict(accession='MANUAL', taxid='999',
        bioproject='P9', reason='misidentified', source_build='historical', reviewer='curator').items()
        if key in fields}
    write_tsv(root / cfg['excluded_accessions'], fields, [manual])
    monkeypatch.setenv('FAKE_GG_LOW_BUSCO', '1')
    finish(build, 1)
    before = (root / cfg['excluded_accessions']).read_bytes()
    preview = catalog.record_successes(root, cfg, build, runs=['A1'], dry_run=True)
    assert preview['excluded_runs'] == ['A1']
    assert (root / cfg['excluded_accessions']).read_bytes() == before
    catalog.record_successes(root, cfg, build, runs=['A1'])
    rows = read_tsv(root / cfg['excluded_accessions'])
    assert all(rows[0][key] == value for key, value in manual.items())
    assert rows[1]['accession'] == 'A1'
    assert rows[1]['source_build'] == 'first'
    assert rows[1].get('reviewer', '') == ''
    assert [row['accession'] for row in rows] == ['MANUAL', 'A1']
    assert read_tsv(root / cfg['accepted_samples']) == []


def test_record_subset_only_excludes_requested_busco_failures(sample_project, monkeypatch):
    root, cfg, build = sample_project
    monkeypatch.setenv('FAKE_GG_LOW_BUSCO', '1')
    for index in (1, 2):
        dataset.worker(build, index)
    result = catalog.record_successes(root, cfg, build, runs=['B1'])
    assert result['excluded_runs'] == ['B1']
    assert [row['accession'] for row in read_tsv(root / cfg['excluded_accessions'])] == ['B1']
    result = catalog.record_successes(root, cfg, build, runs=['A1'])
    assert result['excluded_runs'] == ['A1']
    assert [row['accession'] for row in read_tsv(root / cfg['excluded_accessions'])] == ['B1', 'A1']


def test_busco_passing_sample_waits_for_quantification_before_adoption(sample_project, monkeypatch):
    root, cfg, build = sample_project
    monkeypatch.setenv('FAKE_GG_FAIL_QUANT', '1')
    with pytest.raises(subprocess.CalledProcessError): finish(build, 1)
    result = catalog.record_successes(root, cfg, build, runs=['A1'])
    assert result['recorded_runs'] == result['excluded_runs'] == []
    assert result['samples'][0]['status'] == 'incomplete'
    assert read_tsv(root / cfg['accepted_samples']) == []
    assert read_tsv(root / cfg['excluded_accessions']) == []
    monkeypatch.delenv('FAKE_GG_FAIL_QUANT')
    finish(build, 1)
    result = catalog.record_successes(root, cfg, build, runs=['A1'])
    assert result['recorded_runs'] == ['A1']
    assert result['excluded_runs'] == []


def test_memory_retry_changes_internal_tool_budget_and_preserves_successful_sample(sample_project, monkeypatch):
    root, cfg, build = sample_project
    finish(build, 1)
    monkeypatch.setenv('FAKE_GG_FAIL_ASSEMBLY', '1')
    with pytest.raises(subprocess.CalledProcessError): dataset.worker(build, 2)
    monkeypatch.delenv('FAKE_GG_FAIL_ASSEMBLY')
    resources = root / 'retry.yaml'
    resources.write_text('slurm:\n  per_job_resources:\n    sample:\n      cpus: 8\n      mem_gb: 256\n')
    commands = dataset.submit(build, until='quant', dry_run=True, resources=resources)
    assert '--cpus-per-task=8' in commands[0] and '--mem=256000M' in commands[0]
    assert '--array=2,3' in commands[0]
    monkeypatch.setenv('SLURM_CPUS_PER_TASK', '8')
    monkeypatch.setenv('SLURM_MEM_PER_NODE', '256000')
    monkeypatch.setenv('GG_MEM_TOTAL_GB', '12')
    dataset.worker(build, 2)
    budgets = native_events(build)[-1]['budgets']
    assert budgets == {'GG_TASK_CPUS': '8', 'GG_MEM_TOTAL_GB': '250', 'GG_MEM_TOOL_GB': '246'}
    dataset.worker(build, 1)
    assert sum(e['species'] == 'Alpha_plant' for e in native_events(build)) == 1
    assert read_tsv(root / 'config/excluded_accessions.tsv') == []


def test_damaged_matching_stage_is_reported_before_reuse(sample_project):
    root, cfg, build = sample_project
    dataset.worker(build, 1)
    item = dataset.load(build)['items'][0]
    product = dataset.item_products(dataset.load(build), item)['reference']
    Path(product['cds']['path']).write_bytes(b'damaged')
    config = root / 'config/build.yaml'
    definition = yaml.safe_load(config.read_text()); definition['reuse_from'] = str(build)
    config.write_text(yaml.safe_dump(definition))
    report = dataset.plan(root, config, name='second')[-1]
    assert report[0]['assembly'] == 'conflict'
    assert 'registered file changed' in report[0]['reason']
    with pytest.raises(ValueError, match='resolve conflicts'): dataset.prepare(root, 'second', config)
    assert not (root / 'results/second').exists()


def test_all_reused_partial_build_keeps_container_evidence_and_reviewed_subset(sample_project):
    root, cfg, old = sample_project
    finish(old, 1)
    finish(old, 2)
    rows = read_tsv(root / 'input/metadata.tsv')[:2]
    write_tsv(root / 'input/metadata.tsv', list(rows[0]), rows)
    config = root / 'config/build.yaml'
    definition = yaml.safe_load(config.read_text()); definition['reuse_from'] = 'auto'
    config.write_text(yaml.safe_dump(definition))
    second = dataset.prepare(root, 'second', config)
    assert dataset.submit(second, until='quant', dry_run=True) == []
    result = catalog.record_successes(root, cfg, second, runs=['B1'])
    assert result['recorded_runs'] == ['B1']
    assert [r['run'] for r in read_tsv(root / cfg['accepted_samples'])] == ['B1']
    with pytest.raises(ValueError, match='runs absent from build'):
        catalog.record_successes(root, cfg, second, runs=['UNKNOWN'])


@pytest.mark.parametrize('low_busco', [False, True])
def test_adoption_serializes_with_metadata_updates_when_tables_have_different_directories(
        sample_project, monkeypatch, low_busco):
    from dataset_assets import locked
    root, cfg, build = sample_project
    monkeypatch.setenv('FAKE_GG_LOW_BUSCO', '1' if low_busco else '')
    finish(build, 1)
    before = (root / cfg['excluded_accessions']).read_bytes()
    with locked((root / cfg['previous_metadata']).parent / '.metadata.lock'):
        with pytest.raises(BlockingIOError): catalog.record_successes(root, cfg, build)
    assert not (root / cfg['accepted_samples']).exists()
    assert (root / cfg['excluded_accessions']).read_bytes() == before


def test_changed_assembly_settings_prevent_automatic_reuse(sample_project):
    root, cfg, old = sample_project
    finish(old, 1)
    config = root / 'config/build.yaml'
    definition = yaml.safe_load(config.read_text())
    definition.update(reuse_from='auto')
    definition['genegalleon'].setdefault('settings', {})['assembly_method'] = 'Trinity'
    config.write_text(yaml.safe_dump(definition))
    assert dataset.plan(root, config, name='changed')[-1][0]['assembly'] == 'pending'


def test_auto_reuse_merges_matching_stages_and_records_each_source(sample_project, monkeypatch):
    root, cfg, old = sample_project
    monkeypatch.setenv('FAKE_GG_FAIL_BUSCO', '1')
    with pytest.raises(subprocess.CalledProcessError): dataset.worker(old, 1)
    monkeypatch.delenv('FAKE_GG_FAIL_BUSCO')
    config = root / 'config/build.yaml'
    definition = yaml.safe_load(config.read_text()); definition['reuse_from'] = str(old)
    config.write_text(yaml.safe_dump(definition))
    providers = {}
    for stage in ('busco', 'quant'):
        provider = dataset.prepare(root, stage + '_provider', config)
        dataset.submit(provider, until='quant', dry_run=True)
        if stage == 'busco':
            monkeypatch.setenv('FAKE_GG_FAIL_QUANT', '1')
            with pytest.raises(subprocess.CalledProcessError): dataset.worker(provider, 1)
            monkeypatch.delenv('FAKE_GG_FAIL_QUANT')
        else:
            dataset.worker(provider, 1)
        providers[stage] = provider
    definition['reuse_from'] = 'auto'
    config.write_text(yaml.safe_dump(definition))
    merged = dataset.prepare(root, 'merged', config)
    manifest = dataset.load(merged)
    products = dataset.item_products(manifest, manifest['items'][0])
    assert all(products[key] for key in ('reference', 'busco', 'quant'))
    for stage in ('busco', 'quant'):
        assert products[stage]['provenance']['reused_from']['path'] == str(providers[stage] / 'build.json')


def test_partial_stage_reuse_can_publish_a_completed_database(sample_project, fake_odb,
                                                            frozen_reference, command_environment):
    from build_products import complete, load_complete
    from test_portable_build import execute
    if not shutil.which('snakemake') or not shutil.which('seqkit'):
        pytest.skip('Snakemake and seqkit required')
    root, cfg, old = sample_project
    finish(old, 1)
    finish(old, 2)
    catalog.record_successes(root, cfg, old)
    write_tsv(root / 'config/excluded_accessions.tsv', ['accession', 'reason'],
              [{'accession': 'G1', 'reason': 'unusable_reads'}])
    catalog.update_metadata(root, cfg, 'work/healthy', root / 'input/metadata.tsv')
    config = root / 'config/build.yaml'
    definition = yaml.safe_load(config.read_text()); definition['reuse_from'] = 'auto'
    config.write_text(yaml.safe_dump(definition))
    taxonomy = root / 'resources/taxonomy/taxa.sqlite'; taxonomy.parent.mkdir(parents=True)
    shutil.copy2(root / 'input/taxa.sqlite', taxonomy)
    reference = root / 'resources/orthodb/v12_3193'; reference.parent.mkdir(parents=True)
    reference.symlink_to(frozen_reference, target_is_directory=True)
    env = command_environment({'python': os.sys.executable, 'seqkit': shutil.which('seqkit'), 'ODB-mapper': fake_odb})
    second = dataset.prepare(root, 'healthy', config)
    assert dataset.submit(second, until='quant', dry_run=True) == []
    dataset.materialize(second)
    execute(root, second, 'database', env, mapping=True)
    complete(second)
    assert set(load_complete(second)['products']) == {'Alpha_plant_A1', 'Beta_sp-X_B1'}
    assert catalog.record_successes(root, cfg, second)['recorded_runs'] == ['A1', 'B1']
