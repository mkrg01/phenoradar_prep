#!/usr/bin/env python3
"""Stream original per-sequence TPM for selected samples into one long table."""
import argparse
import csv
import json
import math

from common import atomic_writer, read_tsv, species_from_gene_id, text_open
from sample_identity import sample_run_map


FIELDS = ["species", "sequence_id", "tpm"]


def merge(samples, output):
    rows = sorted(read_tsv(samples), key=lambda row: (row["species"], row["run"]))
    if not rows:
        raise ValueError("sequence expression requires selected samples")
    sample_run_map(rows)
    count = 0
    with atomic_writer(output) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(FIELDS)
        for sample in rows:
            seen = set()
            path = sample["abundance"]
            with text_open(path, newline="") as source:
                reader = csv.DictReader(source, delimiter="\t")
                fields = reader.fieldnames or []
                if len(fields) != len(set(fields)) or not {"target_id", "tpm"} <= set(fields):
                    raise ValueError(f"abundance table requires unique target_id and tpm columns: {path}")
                for row in reader:
                    if None in row or any(value is None for value in row.values()):
                        raise ValueError(f"malformed abundance row: {path}:{reader.line_num}")
                    gene = row["target_id"]
                    if gene in seen:
                        raise ValueError(f"duplicate target_id in {path}: {gene}")
                    if species_from_gene_id(gene) != sample["species"]:
                        raise ValueError(f"target_id belongs to a different sample: {gene}")
                    seen.add(gene)
                    try:
                        value = float(row["tpm"])
                    except ValueError as error:
                        raise ValueError(f"invalid TPM in {path}: {gene}") from error
                    if not math.isfinite(value) or value < 0:
                        raise ValueError(f"TPM must be finite and nonnegative in {path}: {gene}")
                    # Keep the original numeric string, zeros, and every quantified
                    # copy, including targets outside the selected/mapped OG set.
                    writer.writerow([sample["species"], gene, row["tpm"]])
            if not seen:
                raise ValueError(f"empty abundance table: {path}")
            count += len(seen)
    return {"samples": len(rows), "sequences": count}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", required=True)
    parser.add_argument("--output", required=True)
    print(json.dumps(merge(**vars(parser.parse_args()))), flush=True)
