"""Explicit database reuse replaces shared product and ODB cache settings."""
import json
import shutil
from pathlib import Path

import pytest
import yaml

import dataset
from common import read_tsv, write_json, write_tsv
from dataset_assets import digest, identities, record
from database_fixtures import database_from_stages
from test_datasets import dataset_project, imported
from test_portable_build import completed_project, execute
from build_products import complete, load_complete


def configure(root, sources, **values):
    path = root/'config/build.yaml'
    cfg = yaml.safe_load(path.read_text())
    cfg.update(reuse_from=sources, **values)
    path.write_text(yaml.safe_dump(cfg))
    return path


@pytest.mark.parametrize('value', ['', ' ', False, 1, {}, [None], ['valid', 1]])
def test_invalid_reuse_sources(dataset_project, value):
    root = dataset_project
    with pytest.raises(ValueError, match='reuse_from must be'):
        dataset.plan(root, configure(root, value))
    assert not (root/'results').exists()


@pytest.mark.parametrize('key', ['store', 'odb.cache_dir'])
def test_separate_cache_settings_are_rejected(dataset_project, key):
    root = dataset_project
    config = configure(root, None)
    cfg = yaml.safe_load(config.read_text())
    if key == 'store': cfg[key] = 'old/store'
    else: cfg['odb']['cache_dir'] = 'old/odb'
    config.write_text(yaml.safe_dump(cfg))
    with pytest.raises(ValueError, match='unknown build'):
        dataset.plan(root, config)


def test_scalar_list_and_duplicate_paths_have_same_plan_without_writes(dataset_project):
    root = dataset_project
    source = imported(root)
    before = {p: p.read_bytes() for p in source.rglob('*') if p.is_file()}
    reports = []
    for sources in (str(source.relative_to(root)), [str(source)], [str(source), str(source.parent)]):
        cfg, _, _, _, _, _, report = dataset.plan(root, configure(root, sources))
        assert len(cfg['reuse_from']) == 1
        reports.append(report)
    assert reports[0] == reports[1] == reports[2]
    assert {r['mapping'] for r in reports[0]} == {'reuse'}
    assert not (root/'results').exists()
    assert {p: p.read_bytes() for p in before} == before


def test_multiple_databases_supply_union_and_missing_samples_stay_pending(dataset_project):
    root = dataset_project
    imported(root)
    _, items = identities(root/'input/metadata.tsv')
    store = root/'resources/dataset_assets'
    first = database_from_stages(root, store, items[:1], 'first_source')
    second = database_from_stages(root, store, items[1:2], 'second_source')
    config = configure(root, [str(first), str(second)])
    report = dataset.plan(root, config)[-1]
    assert [r['assembly'] for r in report] == ['reuse', 'reuse', 'pending']
    # Missing samples must not cause software acquisition during planning.
    assert not (root/'resources/software').exists()


@pytest.mark.parametrize('reverse', [False, True])
def test_conflicting_duplicate_abundance_is_not_resolved_by_list_order(dataset_project, reverse):
    root = dataset_project
    first = imported(root)
    second = root/'resources/other/database'
    shutil.copytree(first, second)
    path = second/'manifest.json'
    data = json.loads(path.read_text())
    product = data['products']['Alpha_plant_A1']
    relative = product['abundance']['path']
    abundance = second/relative
    rows = read_tsv(abundance); rows[0]['tpm'] = '99'
    write_tsv(abundance, list(rows[0]), rows)
    entry = dict(record(abundance), path=relative)
    product['abundance'] = entry
    data['files'] = [entry if p['path'] == relative else p for p in data['files']]
    data['sha256'] = digest({k:v for k,v in data.items() if k != 'sha256'})
    write_json(path, data)
    sources = [str(first), str(second)]
    if reverse: sources.reverse()
    config = configure(root, sources)
    report = dataset.plan(root, config)[-1]
    assert report[0]['assembly'] == 'conflict'
    assert 'conflicting reuse_from databases for Alpha_plant_A1' in report[0]['reason']
    with pytest.raises(ValueError, match='resolve conflicts'):
        dataset.prepare(root, 'conflicting', config)
    assert not (root/'results/conflicting').exists()


def test_incomplete_source_fails_before_preparation(dataset_project):
    root = dataset_project
    source = root/'results/incomplete/database'; source.mkdir(parents=True)
    with pytest.raises(ValueError, match='build is incomplete'):
        dataset.prepare(root, 'new', configure(root, str(source)))
    assert not (root/'results/new').exists()


def test_prepared_build_has_its_own_cache_and_no_source_dependency(dataset_project):
    root = dataset_project
    first = imported(root)
    second = root/'resources/identical/database'; shutil.copytree(first, second)
    config = configure(root, [str(first), str(second)])
    build = dataset.prepare(root, 'reused', config)
    manifest = dataset.load(build)
    cache = build/'work/cache'
    assert manifest['config']['store'] == str(cache/'products')
    assert manifest['analysis']['odb']['cache_dir'] == str(cache/'odb')
    assert manifest['analysis']['translation_cache'] == str(cache/'proteins')
    assert manifest['analysis']['expression_cache'] == str(cache/'expression')
    assert len(json.loads((build/'reuse.json').read_text())['sources']) == 2
    first.parent.rename(first.parent.with_name('source_offline'))
    second.parent.rename(second.parent.with_name('other_offline'))
    assert dataset.submit(build, until='quant', dry_run=True) == []
    assert len(read_tsv(dataset.materialize(build)/'metadata.tsv')) == 3
    assert not (root/'migrations').exists()


def test_database_reuse_retains_expression_without_mapper_or_translation(completed_project):
    root, source_build, env, events = completed_project
    source = source_build/'database'
    other = root/'results/identical/database'; shutil.copytree(source, other)
    config = configure(root, [str(root/'resources/source/products'), str(source), str(other)])
    first_expression = (source/'expression/runs/A1.tsv').read_bytes()
    build = dataset.prepare(root, 'combined', config)
    assert dataset.submit(build, until='quant', dry_run=True) == []
    # Neither a mapper nor a translator is available to silently redo reused work.
    for name in ('ODB-mapper', 'seqkit'):
        executable = Path(env['PATH'].split(':')[0])/name
        executable.unlink(); executable.write_text('#!/bin/sh\nexit 99\n'); executable.chmod(0o755)
    source_build.rename(source_build.with_name('source_unavailable'))
    other.parent.rename(other.parent.with_name('other_unavailable'))
    dataset.materialize(build)
    execute(root, build, 'database', env, mapping=True)
    complete(build)
    assert (build/'database/expression/runs/A1.tsv').read_bytes() == first_expression
    assert (build/'database/expression/runs/A1.tsv').samefile(root/'results/source_unavailable/database/expression/runs/A1.tsv')
    assert all(json.loads(p.read_text())['reused'] for p in (build/'work/database/proteins').glob('*.json'))
    assert load_complete(build)['products'].keys() == load_complete(root/'results/source_unavailable')['products'].keys()
    assert not events.exists()
    assert (build/'database/provenance/reuse.json').is_file()


def test_two_sources_and_one_missing_sample_only_compute_missing_work(completed_project):
    from test_datasets import fake_genegalleon, native_events
    root, baseline, env, events = completed_project
    frozen = dataset.load(baseline)
    items = frozen['items']
    first = database_from_stages(root, frozen['config']['store'], items[:1], 'alpha_source')
    second = database_from_stages(root, frozen['config']['store'], items[1:2], 'beta_source')
    fake_genegalleon(root)
    config = configure(root, [str(first), str(second)])
    build = dataset.prepare(root, 'mixed', config)
    commands = dataset.submit(build, until='quant', dry_run=True)
    assert len(commands) == 3 and all('--array=3%64' in cmd for cmd in commands)
    for stage in dataset.STAGES: dataset.worker(build, stage, 3)
    assert len(native_events(build)) == 3
    assert {event['species'] for event in native_events(build)} == {'Gamma_plant'}
    dataset.materialize(build)
    execute(root, build, 'database', env, mapping=True)
    complete(build)
    planned = json.loads((build/'work/database/orthogroups/mapping/incremental_plan/plan.json').read_text())
    assert planned['reused_species'] == ['Alpha_plant_A1', 'Beta_sp-X_B1']
    assert planned['mapped_species'] == ['Gamma_plant_G1']
    assert len(events.read_text().splitlines()) == 1
    assert set(load_complete(build)['products']) == {item['species'] for item in items}
