"""Reused checks remain conditional on content identity and current file attributes."""
import os
from pathlib import Path
import shutil
import threading

import pytest

import dataset_assets as assets
from verification_cache import VerificationCache, verification_cache


@pytest.fixture
def moved(tmp_path, monkeypatch):
    original = tmp_path / "original"
    original.write_bytes(b"original contents")
    entry = assets.record(original)
    target = tmp_path / "transferred"
    shutil.copy2(original, target)
    entry = dict(entry, path=str(target))
    calls = []
    hash_file = assets.sha256
    def counted(path):
        calls.append(Path(path))
        return hash_file(path)
    monkeypatch.setattr(assets, "sha256", counted)
    return target, entry, calls


def test_verified_transfer_is_hashed_once_across_commands(tmp_path, moved):
    target, entry, calls = moved
    previous = dict(entry)
    cache = tmp_path / "cache.sqlite"
    with verification_cache(cache):
        assert assets.verify(entry) == target
        assert assets.verify(entry) == target
    with verification_cache(cache):
        assert assets.verify(entry) == target
    assert calls == [target]
    assert entry == previous  # Immutable manifest records retain their original stat values.


def test_cached_check_detects_same_size_edit_with_restored_mtime(tmp_path, moved):
    target, entry, calls = moved
    cache = tmp_path / "cache.sqlite"
    with verification_cache(cache):
        assets.verify(entry)
    previous = target.stat()
    target.write_bytes(b"modified contents")
    os.utime(target, ns=(previous.st_atime_ns, previous.st_mtime_ns))
    with verification_cache(cache), pytest.raises(ValueError, match="registered file changed"):
        assets.verify(entry)
    assert calls == [target, target]


def test_cached_check_detects_replacement_and_missing_files(tmp_path, moved):
    target, entry, calls = moved
    cache = tmp_path / "cache.sqlite"
    with verification_cache(cache):
        assets.verify(entry)
    # Replace via a separate inode; the saved mtime and bytes still agree.
    replacement = target.with_suffix(".new")
    shutil.copy2(target, replacement)
    replacement.replace(target)
    with verification_cache(cache):
        assets.verify(entry)
        target.unlink()
        with pytest.raises(ValueError, match="registered file missing"):
            assets.verify(entry)
    assert calls == [target, target]


def test_cache_is_bound_to_expected_checksum(tmp_path, moved):
    target, entry, calls = moved
    with verification_cache(tmp_path / "cache.sqlite"):
        assets.verify(entry)
        with pytest.raises(ValueError, match="registered file changed"):
            assets.verify(dict(entry, sha256="0" * 64))
    assert calls == [target, target]


@pytest.mark.parametrize("unavailable", ["corrupt", "parent_is_file"])
def test_unavailable_cache_falls_back_to_verification(tmp_path, moved, unavailable):
    target, entry, calls = moved
    cache = tmp_path / "cache.sqlite"
    if unavailable == "corrupt":
        cache.write_bytes(b"not a SQLite database")
    else:
        parent = tmp_path / "blocked"
        parent.write_bytes(b"not a directory")
        cache = parent / "cache.sqlite"
    with verification_cache(cache):
        assets.verify(entry)
        assets.verify(entry)
    assert calls == [target]


def test_changes_during_hashing_are_not_cached(tmp_path, moved, monkeypatch):
    target, entry, _ = moved
    hash_file = assets.sha256
    def changing(path):
        target.write_bytes(b"modified contents")
        return entry["sha256"]
    monkeypatch.setattr(assets, "sha256", changing)
    cache_path = tmp_path / "cache.sqlite"
    with verification_cache(cache_path) as cache:
        assert assets.verify(entry) == target  # Preserve the original checksum-based result.
        assert not cache.matches(target, entry["sha256"], assets.stat_identity(target))
    monkeypatch.setattr(assets, "sha256", hash_file)
    with verification_cache(cache_path), pytest.raises(ValueError, match="registered file changed"):
        assets.verify(entry)


def test_independent_cache_writers_preserve_each_others_checks(tmp_path):
    path = tmp_path / "cache.sqlite"
    first, second = VerificationCache(path), VerificationCache(path)
    first.remember(tmp_path / "first", "a" * 64, [1, 2, 3, 4, 5])
    second.remember(tmp_path / "second", "b" * 64, [6, 7, 8, 9, 10])
    first.close()
    second.close()
    loaded = VerificationCache(path)
    try:
        assert loaded.matches(tmp_path / "first", "a" * 64, [1, 2, 3, 4, 5])
        assert loaded.matches(tmp_path / "second", "b" * 64, [6, 7, 8, 9, 10])
    finally:
        loaded.close()


def test_bulk_verification_hashes_four_files_concurrently_and_reuses_results(tmp_path, monkeypatch):
    entries = []
    for index in range(4):
        path = tmp_path / str(index)
        path.write_bytes(f'contents {index}'.encode())
        entries.append(dict(assets.record(path), stat=[]))
    barrier = threading.Barrier(4)
    calls = []
    hash_file = assets.sha256
    def simultaneous(path):
        calls.append(Path(path))
        barrier.wait(timeout=10)
        return hash_file(path)
    monkeypatch.delenv('SLURM_CPUS_PER_TASK', raising=False)
    monkeypatch.setattr(assets, 'sha256', simultaneous)
    with verification_cache(tmp_path / 'cache.sqlite'):
        assert assets.verify_many(entries) == [None] * 4
        assert assets.verify_many(entries) == [None] * 4
    assert set(calls) == {Path(entry['path']) for entry in entries}
    assert len(calls) == 4


def test_bulk_verification_shares_aliases_but_detects_ctime_changes(tmp_path, monkeypatch):
    path = tmp_path / 'original'
    path.write_bytes(b'original contents')
    entry = assets.record(path)
    aliases = [path, tmp_path / 'alias1', tmp_path / 'alias2']
    for alias in aliases[1:]:
        os.link(path, alias)
    entries = [dict(entry, path=str(alias), stat=[]) for alias in aliases]
    calls, hash_file = [], assets.sha256
    def counted(path):
        calls.append(Path(path))
        return hash_file(path)
    monkeypatch.setattr(assets, 'sha256', counted)
    with verification_cache(tmp_path / 'cache.sqlite'):
        assert assets.verify_many(entries) == [None] * 3
        assert assets.verify_many(entries) == [None] * 3
        previous = path.stat()
        path.write_bytes(b'modified contents')
        os.utime(path, ns=(previous.st_atime_ns, previous.st_mtime_ns))
        assert all(isinstance(error, assets.ProductIntegrityError) for error in assets.verify_many(entries))
    assert len(calls) == 2


def test_imported_hardlink_proof_survives_new_command_and_still_checks_changes(tmp_path, monkeypatch):
    source, target = tmp_path / 'original', tmp_path / 'staged'
    source.write_bytes(b'original contents')
    entry = assets.record(source)
    cache_path = tmp_path / 'cache.sqlite'
    with verification_cache(cache_path):
        staged = assets.stage_file(entry, target)
    calls, hash_file = [], assets.sha256
    def counted(path):
        calls.append(Path(path))
        return hash_file(path)
    monkeypatch.setattr(assets, 'sha256', counted)
    with verification_cache(cache_path):
        assets.verify(entry)
        assets.verify(staged)
    assert calls == []
    previous = target.stat()
    target.write_bytes(b'modified contents')
    os.utime(target, ns=(previous.st_atime_ns, previous.st_mtime_ns))
    with verification_cache(cache_path), pytest.raises(assets.ProductIntegrityError):
        assets.verify(entry)
    assert calls == [source]


def test_cache_only_verification_rejects_reads_and_slurm_allocation_bounds_workers(tmp_path, moved, monkeypatch):
    _, entry, calls = moved
    with verification_cache(tmp_path / 'cache.sqlite'), assets.cached_verification():
        assert not assets.verification_cached(entry)
        with pytest.raises(assets.UncachedVerification):
            assets.verify(entry)
        with pytest.raises(assets.UncachedVerification):
            assets.verify_many([entry])
    assert calls == []
    monkeypatch.setenv('SLURM_CPUS_PER_TASK', '2')
    assert assets.verification_workers() == 2
    assert assets.verification_workers(1) == 1


def test_staging_detects_changes_during_copy_verification(tmp_path, monkeypatch):
    source, target = tmp_path / 'original', tmp_path / 'copy'
    source.write_bytes(b'original contents')
    entry = assets.record(source)
    monkeypatch.setattr(assets, 'link_file', shutil.copy2)
    hash_file = assets.sha256
    def changing(path):
        actual = hash_file(path)
        target.write_bytes(b'modified contents')
        return actual
    monkeypatch.setattr(assets, 'sha256', changing)
    with verification_cache(tmp_path / 'cache.sqlite') as cache:
        with pytest.raises(assets.ProductIntegrityError, match='changed while staging'):
            assets.stage_file(entry, target)
        assert not cache.matches(target, entry['sha256'], assets.stat_identity(target))


def test_document_cache_requires_source_identity_checksum_and_intact_payload(tmp_path):
    cache_path = tmp_path / 'cache.sqlite'
    source = tmp_path / 'source'
    identity = [1, 2, 3, 4, 5]
    data = {'empty': '', 'values': [None, 0, False, '  ']}
    with verification_cache(cache_path) as cache:
        cache.remember_document(source, 'test-v1', identity, 'a' * 64, data)
    with verification_cache(cache_path) as cache:
        assert cache.document(source, 'test-v1', identity, 'a' * 64) == data
        assert cache.document(source, 'test-v1', [1, 2, 3, 4, 6], 'a' * 64) is None
        assert cache.document(source, 'test-v1', identity, 'b' * 64) is None
        assert cache.document(source, 'test-v2', identity, 'a' * 64) is None
        with cache.connection:
            cache.connection.execute("UPDATE documents SET payload = ?", (b'corrupt',))
        assert cache.document(source, 'test-v1', identity, 'a' * 64) is None


@pytest.mark.parametrize('kind', ['sample-metadata-v1', 'completed-manifest-v1'])
def test_document_corruption_falls_back_to_source(tmp_path, kind):
    from dataset_assets import digest, identities, record
    from common import write_json, write_tsv
    from portable_build import manifest_data, _MANIFESTS
    source = tmp_path / 'source'
    if kind == 'sample-metadata-v1':
        write_tsv(source, ['scientific_name', 'run', 'taxid', 'future'],
                  [{'scientific_name': 'Alpha plant', 'run': 'A1', 'taxid': '42', 'future': ''}])
        read = lambda: identities(source, source=record(source))
    else:
        data = {'kind': 'completed_build', 'schema_version': 5, 'products': {}}
        data['sha256'] = digest(data)
        write_json(source, data)
        read = lambda: manifest_data(source, consume=True, validate=True)
    with verification_cache(tmp_path / 'cache.sqlite') as cache:
        expected = read()
        with cache.connection:
            cache.connection.execute('UPDATE documents SET payload = ?', (b'corrupt',))
        _MANIFESTS.clear()
        assert read() == expected
