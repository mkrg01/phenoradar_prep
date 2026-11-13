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


def test_plan_reads_each_database_manifest_once(dataset_project, monkeypatch):
    import portable_build
    root = dataset_project
    source = imported(root)
    config = configure(root, str(source))
    manifest = source / "manifest.json"
    portable_build._MANIFESTS.clear()
    calls = []
    read_text = Path.read_text
    def counted(path, *args, **kwargs):
        if path == manifest:
            calls.append(path)
        return read_text(path, *args, **kwargs)
    monkeypatch.setattr(Path, "read_text", counted)
    report = dataset.plan(root, config)[-1]
    assert all(row["mapping"] == "reuse" for row in report)
    assert calls == [manifest]


def test_cli_reuses_verified_files_across_plan_commands(dataset_project, monkeypatch, capsys):
    import sys
    import dataset_assets
    root = dataset_project
    source = imported(root)
    moved = root / "resources/transferred/database"
    shutil.copytree(source, moved)
    config = configure(root, str(moved))
    calls = []
    hash_file = dataset_assets.sha256
    def counted(path):
        calls.append(Path(path))
        return hash_file(path)
    monkeypatch.setattr(dataset_assets, "sha256", counted)
    monkeypatch.setattr(sys, "argv", ["dataset.py", "plan", "--root", str(root), "--config", str(config)])
    assert dataset.main() == 0
    first = capsys.readouterr().out
    assert calls
    calls.clear()
    assert dataset.main() == 0
    assert capsys.readouterr().out == first
    assert calls == []
    assert (root / ".cache/verification.sqlite").is_file()
    monkeypatch.setattr(sys, "argv", ["dataset.py", "submit", "--root", str(root), "--config", str(config),
                                     "--name", "cached_submission", "--until", "quant", "--dry-run"])
    assert dataset.main() == 0
    assert calls == []  # Submission does not verify or stage large products locally.
    capsys.readouterr()
    damaged = next((moved / "cds").glob("*.gz"))
    damaged.write_bytes(damaged.read_bytes() + b"changed")
    monkeypatch.setattr(sys, "argv", ["dataset.py", "plan", "--root", str(root), "--config", str(config)])
    assert dataset.main() == 0
    assert 'rebuild: 1' in capsys.readouterr().out
    attention = read_tsv(root / 'reports/plan/needs_attention.tsv')
    assert 'registered file changed' in attention[0]['reason']


def test_manifest_changes_and_mutable_results_do_not_poison_shared_reads(dataset_project):
    import portable_build
    root = dataset_project
    source = imported(root)
    manifest = source / "manifest.json"
    loaded = load_complete(source, verify_files=False)
    first = next(iter(loaded["products"]))
    original_name = loaded["products"][first]["row"]["scientific_name"]
    loaded["products"][first]["row"]["scientific_name"] = "Changed result"
    assert load_complete(source, verify_files=False)["products"][first]["row"]["scientific_name"] == original_name
    portable_build.completion_path(source)  # Populate the resolver's read cache.
    bad = json.loads(manifest.read_text())
    bad["build_id"] = "modified manifest"
    write_json(manifest, bad)
    with pytest.raises(ValueError, match="build completion record changed"):
        load_complete(source, verify_files=False)


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
    # Keep each database internally consistent so the conflict is between sources.
    from aggregate_tpm import aggregate
    samples = root / 'changed-samples.tsv'
    write_tsv(samples, ['species', 'run', 'abundance'],
              [dict(species='Alpha_plant_A1', run='A1', abundance=str(abundance))])
    aggregate(samples, 'A1', second / data['mapping']['path'],
              second / product['expression']['path'], second / product['expression_qc']['path'])
    for key in ('expression', 'expression_qc'):
        relative = product[key]['path']
        entry = dict(record(second / relative), path=relative)
        product[key] = entry
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
    # This test takes historical sources offline; latest-alias management is
    # exercised separately by the build-versioning tests.
    config = configure(root, [str(root/'resources/source/database'), str(source), str(other)], name_mode='fixed')
    first_expression = (source/'expression/runs/A1.tsv.gz').read_bytes()
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
    assert (build/'database/expression/runs/A1.tsv.gz').read_bytes() == first_expression
    assert (build/'database/expression/runs/A1.tsv.gz').samefile(root/'results/source_unavailable/database/expression/runs/A1.tsv.gz')
    assert all(json.loads(p.read_text())['reused'] for p in (build/'work/database/proteins').glob('*.json'))
    assert load_complete(build)['products'].keys() == load_complete(root/'results/source_unavailable')['products'].keys()
    assert not events.exists()
    assert (build/'database/provenance/reuse.json').is_file()


def test_two_sources_and_one_missing_sample_only_compute_missing_work(completed_project):
    from test_datasets import fake_genegalleon, native_events
    from common import sha256
    root, baseline, env, events = completed_project
    frozen = dataset.load(baseline)
    items = frozen['items']
    reference = sha256(root / 'resources/orthodb/v12_3193/reference.json')
    first = database_from_stages(root, frozen['config']['store'], items[:1], 'alpha_source', reference)
    second = database_from_stages(root, frozen['config']['store'], items[1:2], 'beta_source', reference)
    fake_genegalleon(root)
    config = configure(root, [str(first), str(second)])
    build = dataset.prepare(root, 'mixed', config)
    commands = dataset.submit(build, until='quant', dry_run=True)
    assert len(commands) == 1 and '--array=3' in commands[0]
    dataset.worker(build, 3)
    assert len(native_events(build)) == 1
    assert {event['species'] for event in native_events(build)} == {'Gamma_plant'}
    dataset.materialize(build)
    execute(root, build, 'database', env, mapping=True)
    complete(build)
    planned = json.loads((build/'work/database/orthogroups/mapping/incremental_plan/plan.json').read_text())
    assert planned['reused_species'] == ['Alpha_plant_A1', 'Beta_sp-X_B1']
    assert planned['mapped_species'] == ['Gamma_plant_G1']
    assert len(events.read_text().splitlines()) == 1
    assert set(load_complete(build)['products']) == {item['species'] for item in items}


def test_same_run_with_changed_name_requires_review(dataset_project):
    """A renamed sample must not silently look like an unprocessed accession."""
    from database_reuse import select
    root = dataset_project
    source = imported(root)
    cfg, _ = dataset.settings(root, configure(root, str(source)))
    rows = read_tsv(root / 'input/metadata.tsv')
    rows[0]['scientific_name'] = 'Renamed plant'
    for key in ('species_id', 'analysis_sample_id'):
        rows[0].pop(key, None)
    metadata = root / 'renamed.tsv'
    write_tsv(metadata, list(rows[0]), rows)
    _, items = identities(metadata)
    selected, errors, _ = select(cfg['reuse_from'], items, cfg)
    renamed = next(i for i in items if i['row']['scientific_name'] == 'Renamed plant')
    assert renamed['species'] not in selected
    assert 'changed sample identity requires a frozen taxonomy snapshot' in errors[renamed['species']]
    assert renamed['row']['run'] in errors[renamed['species']]
    report = dataset.plan(root, root / 'config/build.yaml', metadata=metadata, preview_taxonomy=False)[-1]
    attention = next(row for row in report if row['species'] == renamed['species'])
    assert attention['assembly'] == 'conflict'
    assert 'changed sample identity requires a frozen taxonomy snapshot' in attention['reason']

@pytest.mark.parametrize("merged", [False, True])
@pytest.mark.parametrize("defer", [False, True])
def test_taxonomy_change_stages_relabelled_upstream_products(dataset_project, merged, defer):
    import gzip
    import sqlite3
    from build_dispatch import verify_and_stage

    root = dataset_project
    source = imported(root)
    before = load_complete(source)['products']['Alpha_plant_A1']
    taxonomy = root / 'input/taxa.sqlite'
    rows = read_tsv(root / 'input/metadata.tsv')
    old_taxid = rows[0]['taxid']
    new_taxid = '9999' if merged else old_taxid
    with sqlite3.connect(taxonomy) as db:
        if merged:
            db.execute('UPDATE species SET taxid=? WHERE taxid=?', (new_taxid, old_taxid))
            db.execute('INSERT INTO merged VALUES (?, ?)', (old_taxid, new_taxid))
        db.execute("UPDATE species SET spname='Renamed plant' WHERE taxid=?", (new_taxid,))
    config = configure(root, str(source), taxonomy={'database': str(taxonomy), 'cache_dir': None,
                                                    'missing_taxid': 'error'})
    rows[0].update(scientific_name='Renamed plant', taxid=new_taxid)
    metadata = root / 'renamed.tsv'
    write_tsv(metadata, list(rows[0]), rows)
    report = dataset.plan(root, config, metadata=metadata)[-1]
    renamed, = [row for row in report if row['run'] == 'A1']
    assert (renamed['assembly'], renamed['busco'], renamed['quant'], renamed['mapping']) == (
        'reuse', 'reuse', 'reuse', 'pending_inputs')
    build = dataset.prepare(root, 'renamed', config, metadata=metadata, defer=defer)
    manifest = dataset.load(build)
    if defer:
        verify_and_stage(build, manifest)
        manifest = dataset.load(build)
    item, = [item for item in manifest['items'] if item['row']['run'] == 'A1']
    products = dataset.item_products(manifest, item)
    assert item['reference_id'] == products['reference']['reference_id']
    assert products['reference']['taxid'] == new_taxid
    with gzip.open(products['reference']['cds']['path'], 'rt') as current, gzip.open(before['cds']['path'], 'rt') as old:
        assert current.read().replace('Renamed_plant_A1_g', 'Alpha_plant_A1_g') == old.read()
    with gzip.open(products['quant']['abundance']['path'], 'rt') as current, gzip.open(before['abundance']['path'], 'rt') as old:
        assert current.read().replace('Renamed_plant_A1_g', 'Alpha_plant_A1_g') == old.read()
    assert products['busco']['counts'] == before['counts']
    assert products['reference']['provenance']['taxonomy_relabel']['old']['taxid'] == old_taxid
    assert record(before['cds']['path'])['sha256'] == before['cds']['sha256']


@pytest.mark.parametrize('confirmed', [True, False])
def test_taxonomy_reassignment_rebuilds_only_confirmed_run(dataset_project, monkeypatch, confirmed):
    root = dataset_project
    source = imported(root)
    cfg = configure(root, str(source), taxonomy={'database': str(root / 'input/taxa.sqlite'),
                                                  'cache_dir': None, 'missing_taxid': 'error'})
    rows = read_tsv(root / 'input/metadata.tsv')
    rows[0].update(scientific_name=rows[1]['scientific_name'], taxid=rows[1]['taxid'])

    def confirm(self, run, species, taxonomy):
        assert run == 'A1' and species == rows[1]['taxid']
        if not confirmed:
            raise ValueError('BioSample organism does not resolve to the requested species')
        return {'run': run, 'species_taxid': species, 'biosample': 'SAMN12345'}

    monkeypatch.setattr('ncbi_sample.NCBISamples.confirm', confirm)
    metadata = root / 'reidentified.tsv'
    write_tsv(metadata, list(rows[0]), rows)
    report = dataset.plan(root, cfg, metadata=metadata)[-1]
    changed, = [row for row in report if row['run'] == 'A1']
    assert {row['assembly'] for row in report if row['run'] != 'A1'} == {'reuse'}
    if not confirmed:
        assert changed['assembly'] == 'conflict'
        assert 'BioSample organism does not resolve' in changed['reason']
        return
    assert changed['assembly'] == 'pending' and changed['rebuild']
    assert changed['verification'] == 'not_required'
    build = dataset.prepare(root, 'reassigned', cfg, metadata=metadata, defer=True)
    receipt = json.loads((build / 'reuse.json').read_text())
    assert receipt['taxonomy_reassignments']['Beta_sp-X_A1']['biosample'] == 'SAMN12345'
    candidates = json.loads((build / 'reuse_candidates.json').read_text())
    assert set(candidates['rebuilds']) == {'Beta_sp-X_A1'}
    assert len(candidates['samples']) == 2
