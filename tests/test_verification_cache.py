"""Reused checks remain conditional on content identity and current file attributes."""
import os
from pathlib import Path
import shutil

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
