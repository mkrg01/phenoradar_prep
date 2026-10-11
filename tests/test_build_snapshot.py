"""Frozen snapshots preserve metadata while avoiding repeated empty fields."""
import json
import os

import pytest

import dataset
from common import file_record, read_tsv, write_json, write_tsv
from dataset_assets import identities
from test_datasets import dataset_project
from verification_cache import verification_cache


def test_frozen_metadata_roundtrip_and_updated_empty_columns(dataset_project):
    root = dataset_project
    metadata = root / 'input/metadata.tsv'
    rows = read_tsv(metadata)
    for row in rows:
        row.update({f'future_{i}': '' for i in range(200)})
        row['whitespace'] = '  '
    write_tsv(metadata, list(rows[0]), rows)
    with verification_cache(root / '.cache/verification.sqlite'):
        first = dataset.prepare(root, 'sparse', root / 'config/build.yaml', defer=True)
        manifest = dataset.load(first)
        assert [item['row'] for item in manifest['items']] == [item['row'] for item in identities(metadata)[1]]
        raw = json.loads((first / 'build.json').read_text())
        assert raw['schema_version'] == 3
        assert (first / 'build.json').stat().st_size < len(json.dumps(manifest, indent=2).encode())
        assert read_tsv(first / 'metadata.tsv')[0]['future_0'] == ''
        rows[0]['future_0'] = 'now populated'
        for row in rows:
            row['added_column'] = 'new metadata'
        write_tsv(metadata, list(rows[0]), rows)
        second = dataset.prepare(root, 'updated', root / 'config/build.yaml', defer=True)
        updated = dataset.load(second)
        assert updated['items'][0]['row']['future_0'] == 'now populated'
        assert updated['items'][0]['row']['added_column'] == 'new metadata'
        assert dataset.load(first)['items'][0]['row']['future_0'] == ''
        assert 'added_column' not in dataset.load(first)['fields']


def test_legacy_snapshot_index_is_optional_and_detects_restored_mtime_edit(dataset_project):
    root = dataset_project
    build = dataset.prepare(root, 'legacy', root / 'config/build.yaml', defer=True)
    manifest = dataset.load(build)
    manifest['schema_version'] = 2
    filename = build / 'build.json'
    write_json(filename, manifest)
    hashes = json.loads((build / 'checksums.json').read_text())
    hashes['build.json'] = file_record(filename)['sha256']
    write_json(build / 'checksums.json', hashes)
    cache_path = root / '.cache/verification.sqlite'
    with verification_cache(cache_path) as cache:
        assert dataset.load(build) == manifest
        with cache.connection:
            cache.connection.execute("UPDATE documents SET payload = ? WHERE kind = 'prepared-build-v1'", (b'corrupt',))
        assert dataset.load(build) == manifest
    with verification_cache(cache_path):
        mutable = dataset.load(build)
        mutable['items'][0]['row']['run'] = 'changed only in memory'
        assert dataset.load(build) == manifest
        previous = filename.stat()
        filename.write_bytes(filename.read_bytes().replace(b'"schema_version": 2', b'"schema_version": 9', 1))
        os.utime(filename, ns=(previous.st_atime_ns, previous.st_mtime_ns))
        with pytest.raises(ValueError, match='frozen dataset changed: build.json'):
            dataset.load(build)
