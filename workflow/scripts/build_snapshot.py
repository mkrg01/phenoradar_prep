"""Compact frozen build metadata, with lossless support for legacy snapshots."""
from verification_cache import compact_rows, expand_rows


def compact_build(manifest):
    columns, rows = compact_rows(item['row'] for item in manifest['items'])
    return dict(manifest, schema_version=3, row_columns=columns,
                items=[dict(item, row=row) for item, row in zip(manifest['items'], rows)])


def expand_build(manifest):
    columns = manifest.pop('row_columns', None)
    if columns is not None:
        rows = expand_rows(columns, [item['row'] for item in manifest['items']])
        manifest['items'] = [dict(item, row=row) for item, row in zip(manifest['items'], rows)]
    return manifest
