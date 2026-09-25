"""Read and index checkpoint metadata once per file version."""
import csv
from functools import lru_cache
from pathlib import Path


@lru_cache(maxsize=8)
def _read(path, identity):
    with open(path) as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    species, runs = {}, {}
    for row in rows:
        species.setdefault(row['odb_species'], []).append(row)
        runs.setdefault(row['run'], []).append(row)
    return rows, species, runs


def sample_table(path):
    path = Path(path).resolve()
    stat = path.stat()
    return _read(str(path), (stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns))
