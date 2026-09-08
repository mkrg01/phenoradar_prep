#!/usr/bin/env python3
"""Translate CDS with seqkit, preserving and validating the original IDs."""
import argparse
import gzip
import subprocess

from common import atomic_writer, file_record, now, write_json


def fasta_ids(path, compressed=False):
    opener = gzip.open if compressed else open
    ids = []
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                fields = line[1:].split()
                if not fields:
                    raise ValueError(f"empty FASTA identifier: {path}")
                ids.append(fields[0])
    if not ids or len(ids) != len(set(ids)):
        raise ValueError(f"FASTA must contain unique, nonempty identifiers: {path}")
    return ids


def translate(cds, output, provenance, seqkit="seqkit", table=1, threads=1):
    ids = fasta_ids(cds, compressed=str(cds).endswith(".gz"))
    with atomic_writer(output) as handle:
        subprocess.run([seqkit, "translate", "--threads", str(threads), "--transl-table", str(table), cds],
                       stdout=handle, check=True)
    if fasta_ids(output) != ids:
        raise ValueError("translation changed FASTA identifiers")
    version = subprocess.check_output([seqkit, "version"], text=True).strip()
    write_json(provenance, {"created_at": now(), "cds": file_record(cds), "protein": file_record(output),
                            "seqkit": version, "translation_table": table, "sequences": len(ids)})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ["cds", "output", "provenance"]:
        parser.add_argument(f"--{flag}", required=True)
    parser.add_argument("--seqkit", default="seqkit")
    parser.add_argument("--table", type=int, default=1)
    parser.add_argument("--threads", type=int, default=1)
    translate(**vars(parser.parse_args()))
