"""Small shared helpers; no workflow-engine dependency."""
import csv
import gzip
import hashlib
import io
import json
import os
import re
import tempfile
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path


# Fixed across representative/pair selection and phylogenetic inference.
RANDOM_SEED = 12345


def species_from_gene_id(gene_id):
    """Recover the computational sample key from a {sample}_g{number} ID.

    The function name is retained for existing consumers; this is not taxonomy."""
    match = re.fullmatch(r"([A-Za-z0-9][A-Za-z0-9_.-]*)_g[0-9]+", gene_id)
    if match is None:
        raise ValueError(f"gene ID must use {{species}}_g{{number}}: {gene_id!r}")
    return match.group(1)


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def file_record(path):
    path = Path(path).resolve()
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": sha256(path)}


def text_open(path, mode="rt", **kwargs):
    opener = gzip.open if str(path).endswith(".gz") else open
    return opener(path, mode, **kwargs)


def existing_product(path, default_gzip=False):
    """Resolve a compressed product, accepting historical plain products."""
    path = Path(path)
    compressed = Path(str(path) + ".gz")
    return compressed if compressed.is_file() or (default_gzip and not path.is_file()) else path


def content_sha256(entry):
    """Scientific identity is independent of the product's storage encoding."""
    return entry.get("content_sha256", entry["sha256"])


def protein_provenance(path):
    path = Path(path)
    return (path.with_suffix("") if path.suffix == ".gz" else path).with_suffix(".json")


def now():
    return datetime.now(timezone.utc).isoformat()


@contextmanager
def atomic_writer(path, mode="w"):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        if "b" not in mode and path.suffix == ".gz":
            with os.fdopen(fd, "wb") as raw:
                with gzip.GzipFile(filename="", fileobj=raw, mode="wb", compresslevel=6, mtime=0) as zipped:
                    with io.TextIOWrapper(zipped, encoding="utf-8", newline="") as handle:
                        yield handle
        else:
            kwargs = {} if "b" in mode else {"encoding": "utf-8", "newline": ""}
            with os.fdopen(fd, mode, **kwargs) as handle:
                yield handle
        os.replace(name, path)
    finally:
        if os.path.exists(name):
            os.unlink(name)


def write_json(path, value, compact=False):
    with atomic_writer(path) as handle:
        if compact:
            handle.write(json.dumps(value, separators=(",", ":"), ensure_ascii=False, allow_nan=False))
        else:
            json.dump(value, handle, indent=2, ensure_ascii=False, allow_nan=False)
        handle.write("\n")


def read_tsv(path):
    with text_open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, fields, rows):
    with atomic_writer(path) as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
