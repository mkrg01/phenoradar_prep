#!/usr/bin/env python3
"""Verify every reference checksum without downloading or changing any data."""
import argparse
import json
from pathlib import Path

from common import sha256


def verify(reference):
    record = json.loads(Path(reference).read_text())
    inventory = record["inventory"]
    if sha256(inventory["path"]) != inventory["sha256"]:
        raise ValueError("reference inventory checksum mismatch")
    entries = json.loads(Path(inventory["path"]).read_text())
    for entry in entries:
        path = Path(record["data_dir"]) / entry["relative_path"]
        if not path.is_file() or sha256(path) != entry["sha256"]:
            raise ValueError(f"reference checksum mismatch: {path}")
    print(f"Verified {len(entries)} reference files")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", required=True)
    verify(**vars(parser.parse_args()))
