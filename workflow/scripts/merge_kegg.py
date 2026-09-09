#!/usr/bin/env python3
"""Merge selected KO annotations and per-run TPM, retaining missing evidence."""
import argparse
import csv
import json
import math
import re
import sqlite3
import tempfile
from pathlib import Path

from aggregate_ko_tpm import (GENE_FIELDS, HIT_FIELDS, KO_FIELDS, QC_FIELDS, STATUSES,
                              abundance_values, load_annotation, number, selected_samples,
                              summed, table, verify_file)
from common import atomic_writer, write_tsv


def merge(samples, annotation_dir, run_dir, outdir):
    manifest = selected_samples(samples)
    species = {row["species"]: row["odb_species"] for row in manifest}
    output = Path(outdir)
    output.mkdir(parents=True, exist_ok=True)
    sources, all_kos = {}, set()
    reports = []
    # Gene IDs can number in the millions. Enforce global uniqueness on disk.
    with tempfile.TemporaryDirectory(prefix=".kegg-merge-", dir=output) as work:
        with sqlite3.connect(str(Path(work) / "genes.sqlite")) as database:
            database.execute("CREATE TABLE genes (gene_id TEXT PRIMARY KEY, species TEXT, ko TEXT)")
            database.execute("CREATE TABLE support (species TEXT, run TEXT, ko TEXT, tpm_sum TEXT, "
                             "annotated_genes INTEGER, quantified_genes INTEGER, PRIMARY KEY (run, ko))")
            with atomic_writer(output / "genes.tsv") as gene_handle, atomic_writer(output / "gene_kos.tsv") as hit_handle:
                gene_writer = csv.DictWriter(gene_handle, fieldnames=GENE_FIELDS, delimiter="\t", lineterminator="\n")
                hit_writer = csv.DictWriter(hit_handle, fieldnames=HIT_FIELDS, delimiter="\t", lineterminator="\n")
                gene_writer.writeheader()
                hit_writer.writeheader()
                for name, directory in sorted(species.items()):
                    genes, hits, sources[name] = load_annotation(Path(annotation_dir) / directory, name)
                    try:
                        database.executemany("INSERT INTO genes VALUES (?, ?, ?)",
                                             ((gene, name, row["selected_ko"]) for gene, row in genes.items()))
                    except sqlite3.IntegrityError as error:
                        raise ValueError("gene IDs are not globally unique across selected species") from error
                    gene_writer.writerows({key: row[key] for key in GENE_FIELDS} for row in genes.values())
                    hit_writer.writerows({key: row[key] for key in HIT_FIELDS} for row in hits)
            database.execute("CREATE INDEX genes_species ON genes(species)")
            database.execute("CREATE TEMP TABLE targets (gene_id TEXT PRIMARY KEY)")
            for sample in manifest:
                run, name = sample["run"], sample["species"]
                root = Path(run_dir)
                rows = table(root / f"{run}.tsv", KO_FIELDS)
                report = json.loads((root / f"{run}.qc.json").read_text())
                if (report.get("run") != run or report.get("species") != name
                        or not set(QC_FIELDS).issubset(report)):
                    raise ValueError(f"per-run QC does not match manifest/schema: {run}")
                if report["ambiguity"] not in {"drop", "error"}:
                    raise ValueError(f"invalid ambiguity policy in QC: {run}")
                verify_file(root / f"{run}.tsv", report.get("result", {}))
                verify_file(sample["abundance"], report.get("abundance", {}))
                source = sources[name]
                if any(report.get("annotation_provenance", {}).get(key) != source[key]
                       for key in ("sha256", "bytes")):
                    raise ValueError(f"annotation changed after KO aggregation: {run}")
                values = abundance_values(sample["abundance"], run)
                database.execute("DELETE FROM targets")
                database.executemany("INSERT INTO targets VALUES (?)", ((gene,) for gene in values))
                foreign = database.execute("SELECT gene_id FROM targets JOIN genes USING(gene_id) "
                                           "WHERE species != ? LIMIT 1", (name,)).fetchone()
                if foreign:
                    raise ValueError(f"{run}: target ID belongs to a different species: {foreign[0]}")
                expected = dict(database.execute("SELECT ko, count(*) FROM genes WHERE species = ? AND ko != '' GROUP BY ko", (name,)))
                by_ko = {}
                for row in rows:
                    ko = row["ko"]
                    if row["run"] != run or row["species"] != name:
                        raise ValueError(f"per-run table does not match manifest: {run}")
                    if not re.fullmatch(r"K\d{5}", ko) or ko in by_ko:
                        raise ValueError(f"invalid or duplicate KO row in {run}: {ko}")
                    for key in ("annotated_genes", "quantified_genes"):
                        if not re.fullmatch(r"0|[1-9]\d*", row[key]):
                            raise ValueError(f"invalid {key} in {run}: {ko}")
                    annotated, quantified = int(row["annotated_genes"]), int(row["quantified_genes"])
                    if annotated != expected.get(ko) or quantified > annotated:
                        raise ValueError(f"KO support disagrees with annotation: {run}, {ko}")
                    if quantified:
                        number(row["tpm_sum"], f"KO TPM in {run}")
                    elif row["tpm_sum"] != "":
                        raise ValueError(f"unquantified KO must have unavailable TPM: {run}, {ko}")
                    by_ko[ko] = row
                if set(by_ko) != set(expected):
                    raise ValueError(f"per-run KO set disagrees with annotation: {run}")
                numeric = [row for row in rows if int(row["quantified_genes"])]
                retained = summed(number(row["tpm_sum"], "KO TPM") for row in numeric)
                if (report["targets"] != len(values) or report["annotated_kos"] != len(rows)
                        or report["quantified_kos"] != len(numeric)
                        or report["retained_targets"] != sum(int(row["quantified_genes"]) for row in rows)
                        or not math.isclose(number(report["total_tpm"], "QC total TPM"), summed(values.values()))
                        or not math.isclose(number(report["retained_tpm"], "QC retained TPM"), retained)):
                    raise ValueError(f"KO table or input disagrees with QC: {run}")
                partition = summed(number(report[f"{status}_tpm"], "QC partition TPM")
                                   for status in (*STATUSES, "no_protein"))
                if not math.isclose(partition, number(report["total_tpm"], "QC total TPM")):
                    raise ValueError(f"TPM partitions disagree with QC total: {run}")
                database.executemany("INSERT INTO support VALUES (?, ?, ?, ?, ?, ?)",
                                     (tuple(row[key] for key in KO_FIELDS) for row in rows))
                reports.append(report)
                all_kos.update(by_ko)
            write_tsv(output / "ko_support.tsv", KO_FIELDS,
                      (dict(zip(KO_FIELDS, row)) for row in database.execute(
                          "SELECT * FROM support ORDER BY species, run, ko")))
            fields = ["species", "run", "ko", "tpm_sum"]
            write_tsv(output / "ko_tpm_sum.tsv", fields,
                      (dict(zip(fields, row)) for row in database.execute(
                          "SELECT species, run, ko, tpm_sum FROM support WHERE quantified_genes > 0 "
                          "ORDER BY species, run, ko")))
            kos = sorted(all_kos)

            def wide_rows():
                for sample in manifest:
                    run = sample["run"]
                    observed = dict(database.execute("SELECT ko, tpm_sum FROM support WHERE run = ?", (run,)))
                    row = {"species": sample["species"], "run": run}
                    row.update({ko: observed.get(ko, "") for ko in kos})
                    yield row

            write_tsv(output / "ko_tpm_sum_wide.tsv", ["species", "run", *kos], wide_rows())
    write_tsv(output / "mapping_qc.tsv", QC_FIELDS,
              ({key: row[key] for key in QC_FIELDS} for row in reports))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ["samples", "annotation-dir", "run-dir", "outdir"]:
        parser.add_argument(f"--{flag}", required=True)
    merge(**vars(parser.parse_args()))
