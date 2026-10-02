"""Updated snapshots stay identifiable, portable, and separate from retries."""
import json
import os
import shutil
from concurrent.futures import ThreadPoolExecutor, TimeoutError
from threading import Event

import pytest
import yaml

import dataset
from build_products import complete, load_complete
from build_versioning import publish_latest
from common import file_record, read_tsv, write_json, write_tsv
from dataset_assets import digest, locked
from test_datasets import dataset_project, imported
from test_portable_build import execute
from test_submit_cli import invoke, scheduler


def configure(root, provenance=True):
    config = root / 'config/build.yaml'
    cfg = yaml.safe_load(config.read_text())
    cfg.update(name='leaf', name_mode='timestamp')
    if provenance:
        cfg['metadata_provenance'] = 'input/provenance.json'
        write_json(root / cfg['metadata_provenance'], {
            'kind': 'accepted_metadata', 'metadata': file_record(root / 'input/metadata.tsv'),
            'created_at': '2026-10-01T01:00:00+00:00', 'accepted_at': '2026-10-01T02:00:00+00:00',
            'source_stage': {'kind': 'curate', 'source_origin': {'stage': {
                'kind': 'fetch', 'created_at': '2026-09-30T03:00:00+00:00'}}}})
    config.write_text(yaml.safe_dump(cfg))
    return config


def test_timestamp_names_freeze_provenance_and_handle_same_second(dataset_project, monkeypatch):
    root = dataset_project
    imported(root)
    config = configure(root)
    monkeypatch.setattr(dataset, 'now', lambda: '2026-10-02T09:00:00+09:00')
    first = dataset.prepare(root, None, config)
    before = (first / 'build.json').read_bytes()
    frozen_provenance = (first / 'metadata_provenance.json').read_bytes()
    assert first.name == 'leaf_20261002T000000Z'
    # Changed live metadata creates a different snapshot, even within the same second.
    rows = read_tsv(root / 'input/metadata.tsv')[:2]
    write_tsv(root / 'input/metadata.tsv', list(rows[0]), rows)
    configure(root)
    second = dataset.prepare(root, None, config)
    assert second.name == 'leaf_20261002T000000Z_02'
    assert len(dataset.load(first)['items']) == 3
    assert len(dataset.load(second)['items']) == 2
    (root / 'input/provenance.json').unlink()
    saved = dataset.load(first)
    assert saved['dataset_name'] == 'leaf'
    assert saved['prepared_at'] == saved['created_at'] == '2026-10-02T09:00:00+09:00'
    assert saved['metadata_history'] == {'selected_at': '2026-10-01T01:00:00+00:00',
                                        'accepted_at': '2026-10-01T02:00:00+00:00',
                                        'fetched_at': '2026-09-30T03:00:00+00:00'}
    assert saved['config']['metadata_provenance'] == str(first / 'metadata_provenance.json')
    assert (first / 'build.json').read_bytes() == before
    assert (first / 'metadata_provenance.json').read_bytes() == frozen_provenance
    assert not os.path.lexists(root / 'results/leaf_latest')


def test_automatic_submit_creates_new_build_and_path_retry_keeps_snapshot(dataset_project, scheduler, monkeypatch, capsys):
    root = dataset_project
    imported(root)
    config = configure(root, provenance=False)
    monkeypatch.setattr(dataset, 'now', lambda: '2026-10-02T00:00:00+00:00')
    command = ('submit', '--root', root, '--config', config, '--dry-run')
    invoke(monkeypatch, dataset, *command)
    first = root / 'results/leaf_20261002T000000Z'
    before = (first / 'build.json').read_bytes()
    assert f'Resume with --build {first}' in capsys.readouterr().out
    rows = read_tsv(root / 'input/metadata.tsv')[:2]
    write_tsv(root / 'input/metadata.tsv', list(rows[0]), rows)
    invoke(monkeypatch, dataset, *command)
    assert (root / 'results/leaf_20261002T000000Z_02/build.json').is_file()
    (root / 'input/metadata.tsv').unlink()
    config.unlink()
    invoke(monkeypatch, dataset, 'submit', '--root', root, '--build', first, '--dry-run')
    assert (first / 'build.json').read_bytes() == before
    assert scheduler['commands'] == []


def test_explicit_name_keeps_retry_behavior_in_timestamp_mode(dataset_project, scheduler):
    root = dataset_project
    imported(root)
    config = configure(root, provenance=False)
    dataset.submit_named(root, 'reviewed_release', config, dry_run=True)
    build = root / 'results/reviewed_release'
    before = (build / 'build.json').read_bytes()
    with pytest.raises(ValueError, match='reserved'):
        dataset.prepare(root, 'leaf_latest', config)
    (root / 'input/metadata.tsv').unlink()
    dataset.submit_named(root, 'reviewed_release', config, dry_run=True)
    assert (build / 'build.json').read_bytes() == before
    assert dataset.load(build)['dataset_name'] == 'leaf'


@pytest.mark.parametrize('value', ['auto', None, True, 1, {}])
def test_invalid_naming_mode_fails_before_creating_results(dataset_project, value):
    root = dataset_project
    config = configure(root, provenance=False)
    cfg = yaml.safe_load(config.read_text()); cfg['name_mode'] = value
    config.write_text(yaml.safe_dump(cfg))
    with pytest.raises(ValueError, match='name_mode must be'):
        dataset.prepare(root, None, config)
    assert not (root / 'results').exists()


@pytest.mark.parametrize('change', ['stale_metadata', 'missing', 'invalid_path'])
def test_provenance_must_match_metadata(dataset_project, change):
    root = dataset_project
    config = configure(root)
    provenance = root / 'input/provenance.json'
    if change == 'stale_metadata':
        data = json.loads(provenance.read_text()); data['metadata']['sha256'] = '0' * 64
        write_json(provenance, data)
        message, error = 'does not match', ValueError
    elif change == 'missing':
        provenance.unlink()
        message, error = 'provenance.json', FileNotFoundError
    else:
        cfg = yaml.safe_load(config.read_text()); cfg['metadata_provenance'] = True
        config.write_text(yaml.safe_dump(cfg))
        message, error = 'metadata_provenance must be', ValueError
    with pytest.raises(error, match=message):
        dataset.plan(root, config)


@pytest.fixture
def versioned_database(dataset_project, fake_odb, frozen_reference, command_environment, monkeypatch):
    root = dataset_project
    if not shutil.which('snakemake') or not shutil.which('seqkit'):
        pytest.skip('Snakemake and seqkit required')
    imported(root)
    config = configure(root)
    monkeypatch.setattr(dataset, 'now', lambda: '2026-10-02T00:00:00+00:00')
    taxonomy = root / 'resources/taxonomy/taxa.sqlite'; taxonomy.parent.mkdir(parents=True)
    shutil.copy2(root / 'input/taxa.sqlite', taxonomy)
    reference = root / 'resources/orthodb/v12_3193'; reference.parent.mkdir(parents=True)
    reference.symlink_to(frozen_reference, target_is_directory=True)
    env = command_environment({'python': os.sys.executable, 'seqkit': shutil.which('seqkit'), 'ODB-mapper': fake_odb})
    build = dataset.prepare(root, None, config)
    dataset.materialize(build)
    execute(root, build, 'database', env, mapping=True)
    complete(build)
    return root, build


def test_completed_bundle_preserves_dates_and_latest_can_be_repaired(versioned_database):
    root, build = versioned_database
    latest = root / 'results/leaf_latest'
    assert latest.is_symlink() and latest.resolve() == build
    assert os.readlink(latest) == build.name
    data = load_complete(latest)
    assert data['prepared_at'] == '2026-10-02T00:00:00+00:00'
    assert data['created_at'] == data['completed_at']
    assert data['metadata_history']['fetched_at'] == '2026-09-30T03:00:00+00:00'
    source = build / 'database/provenance/metadata_provenance.json'
    assert source.read_bytes() == (root / 'input/provenance.json').read_bytes()
    assert 'provenance/metadata_provenance.json' in {e['path'] for e in json.loads((build / 'database/manifest.json').read_text())['files']}
    latest.unlink()
    complete(build)
    assert latest.resolve() == build
    # Overlapping completions wait for the shared alias rather than failing.
    started = Event()
    def update():
        started.set()
        return publish_latest(build, data)
    with ThreadPoolExecutor(max_workers=1) as executor:
        with locked(root / 'results/.leaf.latest.lock'):
            future = executor.submit(update)
            assert started.wait(timeout=5)
            with pytest.raises(TimeoutError):
                future.result(timeout=0.05)
        assert future.result(timeout=5) == latest
    moved = root.parent / 'copied-database'
    shutil.copytree(build / 'database', moved)
    assert load_complete(moved)['metadata_history'] == data['metadata_history']
    # The latest alias must not duplicate automatic reuse sources.
    cfg = yaml.safe_load((root / 'config/build.yaml').read_text()); cfg['reuse_from'] = 'auto'
    (root / 'config/build.yaml').write_text(yaml.safe_dump(cfg))
    assert dataset.plan(root, root / 'config/build.yaml')[0]['reuse_from'] == [str(build / 'build.json')]


def test_latest_does_not_regress_when_an_older_build_finishes_later(versioned_database, monkeypatch):
    root, build = versioned_database
    monkeypatch.setattr(dataset, 'now', lambda: '2026-10-03T00:00:00+00:00')
    pending = dataset.prepare(root, None, root / 'config/build.yaml')
    with pytest.raises(ValueError, match='mapping incomplete'):
        complete(pending)
    assert (root / 'results/leaf_latest').resolve() == build
    newer = root / 'results/custom_newer'
    shutil.copytree(build / 'database', newer / 'database')
    manifest = newer / 'database/manifest.json'
    data = json.loads(manifest.read_text())
    data.update(build_id=newer.name, prepared_at='2026-10-03T00:00:00+00:00')
    data['sha256'] = digest({k: v for k, v in data.items() if k != 'sha256'})
    write_json(manifest, data)
    publish_latest(newer, load_complete(newer, verify_files=False))
    complete(build)
    assert (root / 'results/leaf_latest').resolve() == newer
    # Preserve an unrelated directory at the reserved alias path.
    latest = root / 'results/leaf_latest'; latest.unlink(); latest.mkdir()
    marker = latest / 'keep.txt'; marker.write_text('keep')
    with pytest.raises(ValueError, match='not a managed build symlink'):
        complete(build)
    assert marker.read_text() == 'keep'
