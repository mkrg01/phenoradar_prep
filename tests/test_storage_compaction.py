"""Compaction retains scientific bytes, native observations, and verified reuse."""
import gzip
import fcntl
import json
from pathlib import Path

import pytest

from common import file_record, write_tsv, write_json
from dataset import worker, submit, load, item_products, materialize, status
from native_archive import compact, verify_archive
from run_odb_chunk import run as run_odb, completed_output, compress_output
from test_storage_cleanup import odb_job
from test_datasets import dataset_project, new_dataset, native_events


def native_tree(root):
    for name, value in {'output/runtime/one.sh': b'run recipe\n' * 100,
                        'output/logs/one.log': b'diagnostic\n' * 100,
                        'output/observations/attempt/run.json': b'{"exit_code":0}',
                        'output/busco/short.txt': b'protected scientific result',
                        'downloads/version_inventory/one.log': b'versions\n' * 100,
                        'input/amalgkit_metadata/sample.tsv': b'original metadata'}.items():
        p = root / name
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_bytes(value)
    (root / 'downloads/empty/subdirectory').mkdir(parents=True)
    return [root / 'output/busco/short.txt']


def test_native_preview_archive_and_repeat_preserve_products_and_observations(tmp_path):
    protected = native_tree(tmp_path)
    before = {str(p.relative_to(tmp_path)): p.read_bytes() for p in tmp_path.rglob('*') if p.is_file()}
    assert compact(tmp_path, protected)['files'] == 3
    assert not (tmp_path / 'diagnostics.tar.gz').exists()
    compact(tmp_path, protected, apply=True)
    saved = json.loads((tmp_path / 'diagnostics.json').read_text())
    verify_archive(tmp_path / 'diagnostics.tar.gz', saved['files'])
    for relative in ('output/observations/attempt/run.json', 'output/busco/short.txt', 'input/amalgkit_metadata/sample.tsv'):
        assert (tmp_path / relative).read_bytes() == before[relative]
    assert not (tmp_path / 'downloads/empty').exists()
    archive = file_record(tmp_path / 'diagnostics.tar.gz')
    assert compact(tmp_path, protected, apply=True)['state'] == 'complete'
    assert file_record(tmp_path / 'diagnostics.tar.gz') == archive


def test_native_interrupted_retirement_retries_after_archive_verification(tmp_path, monkeypatch):
    protected = native_tree(tmp_path)
    unlink = Path.unlink
    victim = tmp_path / 'output/logs/one.log'
    def interrupted(path, *args, **kwargs):
        if path == victim:
            raise PermissionError('interrupted')
        return unlink(path, *args, **kwargs)
    monkeypatch.setattr(Path, 'unlink', interrupted)
    with pytest.raises(PermissionError):
        compact(tmp_path, protected, apply=True)
    assert (tmp_path / 'diagnostics.json').exists() and victim.exists()
    monkeypatch.setattr(Path, 'unlink', unlink)
    compact(tmp_path, protected, apply=True)
    assert not victim.exists()


def test_corrupted_archive_never_retires_remaining_diagnostics(tmp_path):
    protected = native_tree(tmp_path)
    compact(tmp_path, protected, apply=True)
    victim = tmp_path / 'output/logs/one.log'
    victim.parent.mkdir(parents=True)
    victim.write_text('new data')
    (tmp_path / 'diagnostics.tar.gz').write_bytes(b'broken')
    with pytest.raises(ValueError, match='archive changed'):
        compact(tmp_path, protected, apply=True)
    assert victim.read_text() == 'new data'


def test_native_symlinks_are_not_followed(tmp_path):
    protected = native_tree(tmp_path)
    external = tmp_path / 'external.txt'
    external.write_text('caller owned')
    (tmp_path / 'output/link').symlink_to(external)
    compact(tmp_path, protected, apply=True)
    assert external.read_text() == 'caller owned'
    assert (tmp_path / 'output/link').is_symlink()


def test_native_compaction_keeps_normal_worker_reuse_and_materialization(dataset_project):
    build = new_dataset(dataset_project)
    submit(build, until='quant', dry_run=True)
    worker(build, 1)
    manifest = load(build)
    products = item_products(manifest, manifest['items'][0])
    protected = [products['reference']['cds']['path'], products['quant']['abundance']['path'],
                 products['busco']['full']['path'], products['busco']['short']['path']]
    work = build / 'work/genegalleon/New_plant_SRR1'
    diagnostic = work / 'output/runtime/runner.sh'
    diagnostic.parent.mkdir(parents=True, exist_ok=True)
    diagnostic.write_text('native recipe')
    compact(work, protected, apply=True)
    events = native_events(build)
    worker(build, 1)
    assert native_events(build) == events
    assert status(build)[0]['quant'] == 'reuse'
    assert (materialize(build) / 'metadata.tsv').is_file()


def test_verification_retains_native_failure_history(dataset_project):
    from build_dispatch import verify_and_stage
    build = new_dataset(dataset_project)
    submit(build, until='quant', dry_run=True)
    worker(build, 1)
    receipt = build / 'jobs/status/New_plant_SRR1.sample.json'
    original = json.loads(receipt.read_text())
    write_json(receipt, dict(original, state='failed', native_exit_code=2,
                             error='required versions command failed', finished_at='original-native-time'))
    events = native_events(build)
    verify_and_stage(build, load(build))
    saved = json.loads(receipt.read_text())
    assert saved['state'] in {'complete', 'reused'} and saved['verification'] == 'verified'
    assert saved['native_exit_code'] == 2
    assert saved['error'] == 'required versions command failed'
    assert saved['finished_at'] == 'original-native-time'
    assert saved['completion_basis'] == 'verified_products'
    assert native_events(build) == events


@pytest.mark.parametrize('exit_code', [2, None])
def test_compaction_retains_failed_or_legacy_native_work(dataset_project, monkeypatch, exit_code):
    import build_compaction
    build = new_dataset(dataset_project)
    submit(build, until='quant', dry_run=True)
    worker(build, 1)
    receipt = build / 'jobs/status/New_plant_SRR1.sample.json'
    job = json.loads(receipt.read_text())
    job.pop('native_exit_code', None)
    if exit_code is not None:
        job['native_exit_code'] = exit_code
    write_json(receipt, job)
    monkeypatch.setattr(build_compaction, 'load_complete', lambda path: {})
    monkeypatch.setattr(build_compaction, 'jobs', lambda *args: {'sample': {'activity': None}})
    work = build / 'work/genegalleon/New_plant_SRR1'
    diagnostic = work / 'output/logs/failure.log'
    diagnostic.parent.mkdir(parents=True, exist_ok=True)
    diagnostic.write_text('genuine native diagnostic')
    report = build_compaction.run(build, apply=True)
    assert report['samples'][0]['state'] == 'retained'
    assert diagnostic.read_text() == 'genuine native diagnostic'
    assert not (work / 'diagnostics.tar.gz').exists()


def test_compaction_refuses_active_jobs_and_native_lock(dataset_project, monkeypatch):
    import build_compaction
    build = new_dataset(dataset_project)
    submit(build, until='quant', dry_run=True)
    worker(build, 1)
    monkeypatch.setattr(build_compaction, 'load_complete', lambda path: {})
    monkeypatch.setattr(build_compaction, 'jobs', lambda *args: {'sample': {'activity': 'queued_or_running'}})
    with pytest.raises(ValueError, match='active or unresolved'):
        build_compaction.run(build, apply=True)
    monkeypatch.setattr(build_compaction, 'jobs', lambda *args: {'sample': {'activity': None}})
    work = build / 'work/genegalleon/New_plant_SRR1'
    lock_path = work / '.native.lock'
    with lock_path.open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        with pytest.raises(BlockingIOError):
            build_compaction.run(build, apply=True)
    assert not (work / 'diagnostics.tar.gz').exists()


def test_compaction_requires_a_verified_completed_database(dataset_project):
    from build_compaction import run
    build = new_dataset(dataset_project)
    submit(build, until='quant', dry_run=True)
    worker(build, 1)
    with pytest.raises((ValueError, FileNotFoundError)):
        run(build, apply=True)
    assert not list(build.rglob('diagnostics.tar.gz'))


def test_odb_gzip_reuse_corruption_and_legacy_migration(odb_job):
    args, events = odb_job
    original = run_odb(**args)
    expected = {Path(e['path']).name: Path(e['path']).read_bytes() for e in original['results']}
    before = events.read_text()
    compressed = run_odb(**args, gzip_results=True)
    assert events.read_text() == before
    for entry in compressed['results']:
        path = Path(entry['path'])
        value = gzip.decompress(path.read_bytes()) if path.suffix == '.gz' else path.read_bytes()
        assert value == expected[path.name.removesuffix('.gz')]
    assert completed_output(args['output_dir'], original['fingerprint'], args['label']) == compressed
    assert run_odb(**args, gzip_results=True) == compressed
    annotation = args['output_dir'] / 'chunk_000.og.annotations.gz'
    annotation.write_bytes(annotation.read_bytes() + b'corruption')
    assert completed_output(args['output_dir'], original['fingerprint'], args['label']) is None
    run_odb(**args, gzip_results=True)
    assert events.read_text() != before


def test_compressed_chunk_collect_and_incremental_reuse_match_plain(odb_job, tmp_path):
    from incremental_odb import publish, plan
    from mapping_tables import collect, load_tables, read_species
    args, events = odb_job
    args['output_dir'] = tmp_path / 'chunks/chunk_000'
    plain = run_odb(**args)
    samples = tmp_path / 'samples.tsv'
    write_tsv(samples, ['species', 'odb_species'], [{'species': 'Alpha_plant', 'odb_species': 'Alpha_plant'}])
    chunks = tmp_path / 'chunks.json'
    write_json(chunks, [{'chunk': 'chunk_000', 'species': ['Alpha_plant']}])
    proteins = Path(args['manifest']).parent
    plain_output, gzip_output = tmp_path / 'plain_mapping', tmp_path / 'gzip_mapping'
    collect(samples, chunks, args['output_dir'].parent, proteins, plain_output, tmp_path / 'plain_cache')
    compress_output(args['output_dir'], plain['fingerprint'], args['label'])
    collect(samples, chunks, args['output_dir'].parent, proteins, gzip_output, tmp_path / 'gzip_cache')
    collect(samples, chunks, args['output_dir'].parent, proteins, gzip_output, tmp_path / 'gzip_cache')
    assert read_species(plain_output / 'snapshot.json', 'Alpha_plant')[0] == read_species(gzip_output / 'snapshot.json', 'Alpha_plant')[0]
    cache = tmp_path / 'reuse_cache'
    snapshot = publish(samples, args['output_dir'], proteins, cache, args['label'])
    assert (snapshot / 'annotations.tsv.gz').samefile(args['output_dir'] / 'chunk_000.og.annotations.gz')
    assert plan(samples, proteins, tmp_path / 'reuse_plan', cache)['mapped_species'] == []
    assert publish(samples, args['output_dir'], proteins, cache, args['label']) == snapshot
