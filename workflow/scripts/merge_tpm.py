#!/usr/bin/env python3
"""Combine per-run TPM outputs without implicitly pooling biological samples."""
import argparse
import csv
import io
import json
from contextlib import ExitStack
from pathlib import Path

from common import atomic_writer, existing_product, read_tsv


def merge(samples, run_dir, outdir):
    manifest = sorted(read_tsv(samples), key=lambda row: (row["species"], row["run"]))
    fields = ["species", "run", "targets", "protein_genes", "quantified_proteins", "mapped_targets",
              "ambiguous_targets", "retained_targets", "orthogroups", "total_tpm", "mapped_tpm",
              "mapped_tpm_fraction", "retained_tpm", "retained_tpm_fraction", "multimap"]
    values = ("tpm", "tpm_sum")
    names = ("tpm", "tpm_wide", "tpm_sum", "tpm_sum_wide", "mapping_qc")
    with ExitStack() as stack:
        writers = {
            name: csv.writer(stack.enter_context(atomic_writer(Path(outdir) / f"{name}.tsv")),
                             delimiter="\t", lineterminator="\n")
            for name in names
        }
        for value in values:
            writers[value].writerow(["species", "run", "orthogroup", value])
        writers["mapping_qc"].writerow(fields)

        groups, per_run = set(), []
        for sample in manifest:
            run = sample["run"]
            rows = read_tsv(existing_product(Path(run_dir) / f"{run}.tsv"))
            if not rows or any(row["run"] != run or row["species"] != sample["species"] for row in rows):
                raise ValueError(f"per-run table does not match manifest: {run}")
            seen = {row["orthogroup"] for row in rows}
            if len(seen) != len(rows):
                raise ValueError(f"duplicate OG rows in {run}")
            groups.update(seen)
            # Retain only three UTF-8 columns, not a dictionary for every input row.
            # Original numeric strings and CSV quoting survive the round trip.
            with io.StringIO(newline="") as packed:
                packer = csv.writer(packed, delimiter="\t", lineterminator="\n")
                for row in rows:
                    prefix = (row["species"], row["run"], row["orthogroup"])
                    for value in values:
                        writers[value].writerow((*prefix, row[value]))
                    packer.writerow([row["orthogroup"], row["tpm"], row["tpm_sum"]])
                per_run.append(packed.getvalue().encode("utf-8"))
            del rows
            report = json.loads((Path(run_dir) / f"{run}.qc.json").read_text())
            writers["mapping_qc"].writerow([report[key] for key in fields])

        groups = sorted(groups)
        for value in values:
            writers[f"{value}_wide"].writerow(["species", "run", *groups])
        for index, sample in enumerate(manifest):
            packed = per_run[index]
            per_run[index] = None
            with io.StringIO(packed.decode("utf-8"), newline="") as handle:
                by_group = {group: (tpm, summed) for group, tpm, summed in csv.reader(handle, delimiter="\t")}
            del packed
            for column, value in enumerate(values):
                writers[f"{value}_wide"].writerow(
                    [sample["species"], sample["run"],
                     *(by_group.get(group, (0, 0))[column] for group in groups)])
            del by_group


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ["samples", "run-dir", "outdir"]:
        parser.add_argument(f"--{flag}", required=True)
    merge(**vars(parser.parse_args()))
