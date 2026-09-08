#!/usr/bin/env python3
"""Combine per-run TPM outputs without implicitly pooling biological samples."""
import argparse
import json
from pathlib import Path

from common import read_tsv, write_tsv


def merge(samples, run_dir, outdir):
    manifest = sorted(read_tsv(samples), key=lambda row: (row["species"], row["run"]))
    long_rows, reports, per_run = [], [], {}
    groups = set()
    for sample in manifest:
        run = sample["run"]
        rows = read_tsv(Path(run_dir) / f"{run}.tsv")
        if not rows or any(row["run"] != run or row["species"] != sample["species"] for row in rows):
            raise ValueError(f"per-run table does not match manifest: {run}")
        if len({row["orthogroup"] for row in rows}) != len(rows):
            raise ValueError(f"duplicate OG rows in {run}")
        long_rows.extend(rows)
        per_run[run] = {row["orthogroup"]: row for row in rows}
        groups.update(per_run[run])
        reports.append(json.loads((Path(run_dir) / f"{run}.qc.json").read_text()))
    groups = sorted(groups)
    out = Path(outdir)
    for value, name in [("tpm", "tpm"), ("tpm_sum", "tpm_sum")]:
        write_tsv(out / f"{name}.tsv", ["species", "run", "orthogroup", value],
                  ({k: row[k] for k in ["species", "run", "orthogroup", value]} for row in long_rows))
        wide = []
        for sample in manifest:
            run = sample["run"]
            row = {"species": sample["species"], "run": run}
            row.update({og: per_run[run].get(og, {}).get(value, 0) for og in groups})
            wide.append(row)
        write_tsv(out / f"{name}_wide.tsv", ["species", "run", *groups], wide)
    fields = ["species", "run", "targets", "protein_genes", "quantified_proteins", "mapped_targets",
              "ambiguous_targets", "retained_targets", "orthogroups", "total_tpm", "mapped_tpm",
              "mapped_tpm_fraction", "retained_tpm", "retained_tpm_fraction", "multimap"]
    write_tsv(out / "mapping_qc.tsv", fields, ({k: row[k] for k in fields} for row in reports))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ["samples", "run-dir", "outdir"]:
        parser.add_argument(f"--{flag}", required=True)
    merge(**vars(parser.parse_args()))
