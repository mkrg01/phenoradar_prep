#!/usr/bin/env python3
"""Sum OG TPM per run, with explicit ambiguity and renormalization semantics."""
import argparse
import csv
import math
import sqlite3
from collections import defaultdict
from pathlib import Path

from common import file_record, now, read_tsv, write_json, write_tsv


def aggregate(samples, run, database, output, qc, multimap="error"):
    if multimap not in {"error", "drop", "split"}:
        raise ValueError("unknown multimap policy")
    rows = [row for row in read_tsv(samples) if row["run"] == run]
    if len(rows) != 1:
        raise ValueError(f"expected exactly one manifest row for run {run}")
    sample = rows[0]
    values = {}
    with open(sample["abundance"]) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {"target_id", "tpm"}.issubset(reader.fieldnames or []):
            raise ValueError("abundance table requires target_id and tpm")
        for row in reader:
            gene = row["target_id"]
            if not gene or gene in values:
                raise ValueError(f"empty or duplicate target_id in {run}: {gene}")
            value = float(row["tpm"])
            if not math.isfinite(value) or value < 0:
                raise ValueError(f"invalid TPM in {run}: {gene}")
            values[gene] = value
    total = math.fsum(values.values())
    if not values or total <= 0:
        raise ValueError(f"{run}: abundance has no positive TPM")
    db = sqlite3.connect(Path(database).resolve().as_uri() + "?mode=ro", uri=True)
    mapping = defaultdict(list)
    try:
        db.execute("CREATE TEMP TABLE targets (query TEXT PRIMARY KEY)")
        db.executemany("INSERT INTO targets VALUES (?)", ((gene,) for gene in values))
        foreign = db.execute("SELECT genes.query FROM genes JOIN targets USING(query) WHERE species != ? LIMIT 1",
                             (sample["species"],)).fetchone()
        if foreign:
            raise ValueError(f"{run}: target ID belongs to a different species: {foreign[0]}")
        for gene, og in db.execute("SELECT mappings.query, og FROM mappings JOIN targets USING(query)"):
            mapping[gene].append(og)
        protein_genes = db.execute("SELECT count(*) FROM genes WHERE species = ?", (sample["species"],)).fetchone()[0]
        quantified_proteins = db.execute("SELECT count(*) FROM genes JOIN targets USING(query) WHERE species = ?",
                                        (sample["species"],)).fetchone()[0]
    finally:
        db.close()
    ambiguous = [gene for gene, groups in mapping.items() if len(groups) > 1]
    if ambiguous and multimap == "error":
        raise ValueError(f"{run}: {len(ambiguous)} genes map to multiple OGs; choose drop/split explicitly; examples {ambiguous[:5]}")
    summed, included = defaultdict(float), set()
    for gene, groups in mapping.items():
        if len(groups) > 1 and multimap == "drop":
            continue
        for group in groups:
            summed[group] += values[gene] / len(groups)
        included.add(gene)
    matched_tpm = math.fsum(values[gene] for gene in mapping)
    included_tpm = math.fsum(values[gene] for gene in included)
    if included_tpm <= 0:
        raise ValueError(f"{run}: no positive TPM maps to retained OGs; check FASTA/query IDs and ambiguity policy")
    result = [{"species": sample["species"], "run": run, "orthogroup": group,
               "tpm_sum": value, "tpm": value / included_tpm * 1e6}
              for group, value in sorted(summed.items())]
    write_tsv(output, ["species", "run", "orthogroup", "tpm_sum", "tpm"], result)
    write_json(qc, {"created_at": now(), "species": sample["species"], "run": run, "multimap": multimap,
                    "targets": len(values), "protein_genes": protein_genes, "quantified_proteins": quantified_proteins,
                    "mapped_targets": len(mapping), "ambiguous_targets": len(ambiguous), "retained_targets": len(included),
                    "total_tpm": total, "mapped_tpm": matched_tpm, "mapped_tpm_fraction": matched_tpm / total,
                    "retained_tpm": included_tpm, "retained_tpm_fraction": included_tpm / total,
                    "orthogroups": len(result), "abundance": file_record(sample["abundance"])})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ["samples", "run", "database", "output", "qc"]:
        parser.add_argument(f"--{flag}", required=True)
    parser.add_argument("--multimap", choices=["error", "drop", "split"], default="error")
    aggregate(**vars(parser.parse_args()))
