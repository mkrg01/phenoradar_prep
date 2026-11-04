"""Missing taxonomy is excluded before reuse import, arrays, or mapping."""
import json
import sqlite3
from pathlib import Path

import pytest
import yaml

from common import read_tsv
from dataset import load, materialize, plan, prepare, submit
from dataset_assets import ProductIntegrityError
from taxonomy_selection import apply, classify, ensure
from test_datasets import dataset_project, fake_genegalleon, imported


def policy(root, value='exclude'):
    config = root / 'config/build.yaml'
    cfg = yaml.safe_load(config.read_text())
    cfg['taxonomy'].update(database=str(root / 'input/taxa.sqlite'), missing_taxid=value)
    config.write_text(yaml.safe_dump(cfg))
    return config


def remove(root, *taxids):
    with sqlite3.connect(root / 'input/taxa.sqlite') as db:
        db.executemany('DELETE FROM species WHERE taxid=?', [(i,) for i in taxids])


def test_plan_excludes_before_reuse_checks_and_preserves_input(dataset_project, monkeypatch):
    root = dataset_project
    config = policy(root)
    remove(root, 42)
    original = (root / 'input/metadata.tsv').read_bytes()
    inspected = []
    monkeypatch.setattr('dataset.inspect_items', lambda store, items, *a, **kw:
                        inspected.extend(i['row']['run'] for i in items) or [])
    result = plan(root, config)
    assert inspected == ['B1', 'G1']
    assert [i['row']['run'] for i in result[4]] == ['B1', 'G1']
    assert result[5]['excluded'][0]['reason'] == 'taxid_unresolved'
    assert (root / 'input/metadata.tsv').read_bytes() == original
    assert not (root / 'results').exists()


def test_frozen_selection_covers_reused_samples_and_mapping(dataset_project):
    root = dataset_project
    imported(root)
    config = policy(root)
    remove(root, 42)
    build = prepare(root, 'filtered', config, defer=True)
    prepared = (build / 'build.json').read_bytes()
    candidate_metadata = (build / 'metadata.tsv').read_bytes()
    candidates = json.loads((build / 'reuse_candidates.json').read_text())
    omitted = next(v for k, v in candidates['samples'].items() if k.startswith('Alpha_'))
    Path(omitted['product']['cds']['path']).write_bytes(b'excluded damaged product')
    from build_dispatch import verify_and_stage
    verify_and_stage(build, load(build))
    selected = load(build)
    assert [i['row']['run'] for i in selected['items']] == ['B1', 'G1']
    assert selected['excluded'][0]['run'] == 'A1'
    assert selected['excluded'][0]['source'] == 'taxonomy'
    assert (build / 'build.json').read_bytes() == prepared
    assert (build / 'metadata.tsv').read_bytes() == candidate_metadata
    assert not (build / 'jobs/status/Alpha_plant_A1.sample.json').exists()
    assert [r['run'] for r in read_tsv(materialize(build) / 'metadata.tsv')] == ['B1', 'G1']
    jobs = submit(build, dry_run=True)
    assert len(jobs) == 1 and '--array' not in ' '.join(jobs[0])
    data = json.loads((build / 'taxonomy_selection.json').read_text())
    assert data['rows'][0]['action'] == 'exclude'
    assert data['database']['sha256']
    assert Path(selected['analysis']['taxonomy_db']).resolve() == (root / 'input/taxa.sqlite').resolve()


def test_array_indices_and_worker_use_the_same_filtered_order(dataset_project):
    root = dataset_project
    fake_genegalleon(root)
    config = policy(root)
    remove(root, 42)
    build = prepare(root, 'array', config, defer=True)
    from build_dispatch import verify_and_stage
    verify_and_stage(build, load(build))
    commands = submit(build, until='quant', dry_run=True)
    assert '--array=1,2' in commands[0]
    assert [i['row']['run'] for i in load(build)['items']] == ['B1', 'G1']
    from dataset import worker
    worker(build, 1)
    assert (build / 'jobs/status/Beta_sp-X_B1.sample.json').is_file()
    assert not (build / 'jobs/status/Alpha_plant_A1.sample.json').exists()


def test_immediate_prepare_filters_before_product_verification(dataset_project):
    root = dataset_project
    database = imported(root)
    config = policy(root)
    remove(root, 42)
    data = json.loads((database / 'manifest.json').read_text())
    alpha = data['products']['Alpha_plant_A1']['cds']
    (database / alpha['path']).write_bytes(b'damaged excluded sample')
    build = prepare(root, 'immediate', config)
    assert [i['row']['run'] for i in load(build)['items']] == ['B1', 'G1']
    assert not (build / 'jobs/status/Alpha_plant_A1.sample.json').exists()


def test_merged_ids_remain_eligible_without_changing_identity(dataset_project):
    root = dataset_project
    from dataset_assets import identities
    _, items = identities(root / 'input/metadata.tsv')
    with sqlite3.connect(root / 'input/taxa.sqlite') as db:
        db.execute('INSERT INTO merged VALUES (99,42)')
    items[0]['row']['taxid'] = '99'
    rows = classify(items, root / 'input/taxa.sqlite')
    assert rows[0] == dict(species=items[0]['species'], run='A1', taxid='99',
                           resolved_taxid='42', action='include', reason='taxid_merged')
    assert items[0]['row']['taxid'] == '99'


@pytest.mark.parametrize('problem', ['root_missing', 'merged_cycle', 'not_sqlite'])
def test_bad_database_is_an_error_not_an_exclusion(dataset_project, problem):
    root = dataset_project
    from dataset_assets import identities
    _, items = identities(root / 'input/metadata.tsv')
    dbfile = root / 'input/taxa.sqlite'
    if problem == 'not_sqlite':
        dbfile.write_bytes(b'not a database')
    else:
        with sqlite3.connect(dbfile) as db:
            if problem == 'root_missing':
                db.execute('DELETE FROM species WHERE taxid=1')
            else:
                db.execute('DELETE FROM species WHERE taxid=42')
                db.executemany('INSERT INTO merged VALUES (?,?)', [(42,99),(99,42)])
    with pytest.raises((ValueError, sqlite3.DatabaseError)):
        classify(items, dbfile)


@pytest.mark.parametrize('value', ['exclude', 'error'])
def test_all_excluded_or_strict_policy_stops_before_dispatch(dataset_project, value):
    root = dataset_project
    config = policy(root, value)
    fake_genegalleon(root)
    remove(root, 42,43,44)
    build = prepare(root, 'empty', config, defer=True)
    with pytest.raises(ValueError, match='all candidate|taxids absent'):
        ensure(build, load(build))
    assert len(read_tsv(build / 'taxonomy_excluded_runs.tsv')) == 3
    assert not list((build / 'jobs').glob('submission_*.json'))
    with pytest.raises(ValueError, match='all candidate|taxids absent'):
        submit(build, dry_run=True)


def test_changed_snapshot_fails_on_resume_instead_of_reselecting(dataset_project):
    root = dataset_project
    config = policy(root)
    build = prepare(root, 'pinned', config, defer=True)
    ensure(build, load(build))
    remove(root, 42)
    with pytest.raises(ProductIntegrityError):
        load(build)


def test_unverified_candidates_cannot_be_scheduled(dataset_project):
    root = dataset_project
    config = policy(root)
    build = prepare(root, 'pending', config, defer=True)
    with pytest.raises(ValueError, match='taxonomy selection pending'):
        submit(build, dry_run=True)


def test_plan_without_snapshot_is_pending_and_downloads_nothing(dataset_project):
    root = dataset_project
    config = policy(root)
    cfg = yaml.safe_load(config.read_text());cfg['taxonomy']['database']=None
    cfg['taxonomy']['cache_dir']='resources/taxonomy';config.write_text(yaml.safe_dump(cfg))
    result = plan(root, config)
    assert all(r.get('taxonomy') == 'pending' for r in result[-1])
    assert not (root / 'resources/taxonomy').exists()
