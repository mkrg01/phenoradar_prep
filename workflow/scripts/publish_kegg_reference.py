#!/usr/bin/env python3
"""Verify a frozen KEGG snapshot and publish its membership tables for a run."""
import argparse
from pathlib import Path
import shutil

from common import atomic_writer, file_record, now, write_json
from verify_kegg_reference import verify


def publish(reference, qc, modules, pathways):
    root = Path(reference).resolve().parent
    destinations = [Path(path).resolve() for path in [qc, modules, pathways]]
    if len(set(destinations)) != len(destinations):
        raise ValueError("KEGG publication destinations must be distinct")
    if any(path == root or root in path.parents for path in destinations):
        raise ValueError("KEGG publication must not write inside the frozen reference")
    metadata = verify(reference, full=True)
    for source, destination in [(metadata["ko_modules"], modules),
                                (metadata["ko_pathways"], pathways)]:
        with open(source, "rb") as handle, atomic_writer(destination, "wb") as out:
            shutil.copyfileobj(handle, out)
    write_json(qc, {"created_at": now(), "reference_id": metadata["reference_id"],
                    "reference": file_record(reference), "release": metadata["release"],
                    "verification": "all_inventory_checksums",
                    "semantics": "KO membership, not module completeness or activity",
                    "results": [file_record(modules), file_record(pathways)]})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ["reference", "qc", "modules", "pathways"]:
        parser.add_argument(f"--{flag}", required=True)
    publish(**vars(parser.parse_args()))
