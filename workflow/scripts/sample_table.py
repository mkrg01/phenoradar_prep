"""Read and index checkpoint metadata once per file version."""
import csv
from functools import lru_cache
from pathlib import Path

EXPRESSION_BATCH_SIZE = 50


@lru_cache(maxsize=8)
def _read(path, identity):
    with open(path) as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    species, runs = {}, {}
    for row in rows:
        species.setdefault(row.get('odb_species', row.get('species', '')), []).append(row)
        runs.setdefault(row['run'], []).append(row)
    return rows, species, runs


def sample_table(path):
    path = Path(path).resolve()
    stat = path.stat()
    return _read(str(path), (stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns))


@lru_cache(maxsize=8)
def _batches(path, identity, size):
    rows, _, runs = _read(path, identity)
    if any(len(selected) != 1 for selected in runs.values()):
        raise ValueError('expression requires unique runs in the sample manifest')
    rows = sorted(rows, key=lambda row: row['run'])
    batches, by_run = {}, {}
    for offset in range(0, len(rows), size):
        name = f'batch_{offset // size:05d}'
        batches[name] = rows[offset:offset + size]
        by_run.update({row['run']: name for row in batches[name]})
    return batches, by_run


def expression_batches(path, size=EXPRESSION_BATCH_SIZE):
    """Plan stable, bounded batches from the same versioned metadata index."""
    if type(size) is not int or size < 1:
        raise ValueError('expression batch size must be a positive integer')
    path = Path(path).resolve()
    stat = path.stat()
    identity = (stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns)
    return _batches(str(path), identity, size)
