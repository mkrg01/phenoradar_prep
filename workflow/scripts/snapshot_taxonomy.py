#!/usr/bin/env python3
"""Create a local, offline snapshot of an existing ETE taxonomy database."""
import argparse
import os
import sqlite3
import tempfile
from pathlib import Path

from common import file_record, now, write_json


def snapshot(source, destination):
    source, destination = Path(source).resolve(), Path(destination).resolve()
    if not source.is_file():
        raise ValueError(f"source database missing: {source}")
    if destination.exists():
        raise ValueError(f"snapshot already exists; use a new destination to update: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=destination.parent, prefix=".taxonomy.")
    os.close(fd)
    try:
        with sqlite3.connect(source.as_uri() + "?mode=ro", uri=True) as src, sqlite3.connect(tmp) as dst:
            src.backup(dst)
            if dst.execute("PRAGMA quick_check").fetchone()[0] != "ok":
                raise ValueError("invalid taxonomy database")
            dst.execute("SELECT taxid FROM species LIMIT 1")
        os.replace(tmp, destination)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)
    write_json(str(destination) + ".json", {"created_at": now(), "source": str(source),
                                            "snapshot": file_record(destination)})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True)
    parser.add_argument("--destination", required=True)
    snapshot(**vars(parser.parse_args()))
