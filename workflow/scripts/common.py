"""Small shared helpers; no workflow-engine dependency."""
import csv
import hashlib
import json
import os
import tempfile
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def file_record(path):
    path = Path(path).resolve()
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": sha256(path)}


def now():
    return datetime.now(timezone.utc).isoformat()


@contextmanager
def atomic_writer(path, mode="w"):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        kwargs = {} if "b" in mode else {"encoding": "utf-8", "newline": ""}
        with os.fdopen(fd, mode, **kwargs) as handle:
            yield handle
        os.replace(name, path)
    finally:
        if os.path.exists(name):
            os.unlink(name)


def write_json(path, value):
    with atomic_writer(path) as handle:
        json.dump(value, handle, indent=2, ensure_ascii=False, allow_nan=False)
        handle.write("\n")


def read_tsv(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, fields, rows):
    with atomic_writer(path) as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
