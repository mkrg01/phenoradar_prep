#!/usr/bin/env python3
"""Validate input tables, add frozen taxonomy, and select species by BUSCO."""
import argparse
import json
import re
from pathlib import Path

import pandas as pd

from common import atomic_writer, file_record, now, write_json, write_tsv

RANKS = ["kingdom", "phylum", "class", "order", "family", "genus"]
COUNTS = ["busco_cds_single", "busco_cds_duplicated", "busco_cds_fragmented",
          "busco_cds_missing", "busco_cds_total"]


def prepare(metadata, busco, cds_dir, quant_dir, taxonomy_db, outdir,
            threshold=0.5, species_list=None, missing_taxonomy="error"):
    if not 0 <= threshold <= 1:
        raise ValueError("BUSCO threshold must be between 0 and 1")
    meta = pd.read_csv(metadata, sep="\t", dtype=str, keep_default_na=False)
    bus = pd.read_csv(busco, sep="\t", dtype=str, keep_default_na=False)
    for frame, required, label in [(meta, ["scientific_name", "run", "taxid"], "metadata"),
                                   (bus, ["Species", *COUNTS], "BUSCO")]:
        absent = set(required) - set(frame.columns)
        if absent:
            raise ValueError(f"{label}: missing columns: {sorted(absent)}")
        if frame.empty or frame[required].eq("").any().any():
            raise ValueError(f"{label}: empty table or required values")
    if meta["run"].duplicated().any():
        raise ValueError("metadata: run IDs must be unique")
    if bus["Species"].duplicated().any():
        raise ValueError("BUSCO: Species must be unique")
    if (meta.groupby("scientific_name")["taxid"].nunique() > 1).any():
        raise ValueError("metadata: conflicting taxids for the same species")
    reserved = {"species", "odb_species", "busco_percent", *RANKS, *COUNTS, "busco_cds_summary"}
    if reserved.intersection(meta.columns):
        raise ValueError(f"metadata contains derived columns: {sorted(reserved.intersection(meta.columns))}")
    meta["species"] = meta["scientific_name"].str.replace(" ", "_", regex=False)
    meta["odb_species"] = meta["species"].str.replace("-", "_", regex=False)
    for column in ["species", "run"]:
        if not meta[column].map(lambda value: bool(re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", value))).all():
            raise ValueError(f"metadata: unsafe {column} label")
    if meta[["scientific_name", "odb_species"]].drop_duplicates()["odb_species"].duplicated().any():
        raise ValueError("species names collide after space/hyphen normalization")
    for column in COUNTS:
        bus[column] = pd.to_numeric(bus[column], errors="raise")
        if ((bus[column] < 0) | (bus[column] % 1 != 0)).any():
            raise ValueError(f"BUSCO: invalid counts in {column}")
    if (bus["busco_cds_total"] <= 0).any() or not bus[COUNTS[:-1]].sum(axis=1).eq(bus[COUNTS[-1]]).all():
        raise ValueError("BUSCO: counts must sum to a positive total")
    absent = set(meta["scientific_name"]) - set(bus["Species"])
    if absent:
        raise ValueError(f"metadata species missing from BUSCO: {sorted(absent)}")
    joined = meta.merge(bus.rename(columns={"Species": "scientific_name"}),
                        on="scientific_name", how="left", validate="many_to_one")
    joined["busco_percent"] = (joined[COUNTS[0]] + joined[COUNTS[1]]) / joined[COUNTS[-1]] * 100

    # Refuse missing databases BEFORE NCBITaxa can attempt a download or upgrade.
    if not Path(taxonomy_db).is_file():
        raise ValueError(f"taxonomy database missing: {taxonomy_db}; run snapshot_taxonomy.py first")
    from ete4 import NCBITaxa
    ncbi = NCBITaxa(dbfile=str(Path(taxonomy_db).resolve()), update=False)
    taxonomy, unknown = [], []
    try:
        for taxid in sorted(set(joined["taxid"]), key=int):
            try:
                lineage = ncbi.get_lineage(int(taxid))
            except ValueError:
                lineage = []
            if not lineage:
                unknown.append(taxid)
            ranks = ncbi.get_rank(lineage)
            names = ncbi.get_taxid_translator(lineage)
            row = {rank: "" for rank in RANKS}
            for tid in lineage:
                rank = ranks.get(tid)
                if rank in row and not row[rank]:
                    row[rank] = names.get(tid, "")
            taxonomy.append({"taxid": taxid, **row})
    finally:
        ncbi.db.close()
    if unknown and missing_taxonomy == "error":
        raise ValueError(f"taxids absent from frozen taxonomy: {unknown}")
    joined = joined.merge(pd.DataFrame(taxonomy), on="taxid", validate="many_to_one")
    joined["selected"] = (joined[COUNTS[0]] + joined[COUNTS[1]]) / joined[COUNTS[-1]] >= threshold
    requested = None
    if species_list:
        requested = Path(species_list).read_text().splitlines()
        if not requested or any(not value for value in requested) or len(set(requested)) != len(requested):
            raise ValueError("species list must contain unique nonempty species IDs")
        eligible = set(joined.loc[joined["selected"], "species"])
        if set(requested) - eligible:
            raise ValueError(f"requested species absent or below BUSCO threshold: {sorted(set(requested) - eligible)}")
        joined["selected"] &= joined["species"].isin(requested)
    selected = joined.loc[joined["selected"]].sort_values(["species", "run"]).copy()
    if selected.empty:
        raise ValueError("no species passed selection")
    samples = []
    for row in selected.to_dict("records"):
        species, run = row["species"], row["run"]
        cds = (Path(cds_dir) / f"{species}_longestCDS.fa.gz").resolve()
        abundance = (Path(quant_dir) / species / run / f"{run}_abundance.tsv").resolve()
        for path in [cds, abundance]:
            if not path.is_file() or not path.stat().st_size:
                raise ValueError(f"input missing or empty: {path}")
        samples.append({key: row[key] for key in ["scientific_name", "species", "odb_species", "run", "taxid"]}
                       | {"cds": str(cds), "abundance": str(abundance)})
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    for name, frame in [("metadata_all.tsv", joined), ("metadata_high_busco.tsv", selected)]:
        with atomic_writer(out / name) as handle:
            frame.to_csv(handle, sep="\t", index=False)
    write_tsv(out / "samples.tsv", list(samples[0]), samples)
    with atomic_writer(out / "species_high_busco.txt") as handle:
        handle.write("\n".join(sorted(set(selected["species"]))) + "\n")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(joined["busco_percent"], bins=range(102), color="grey")
    ax.axvline(threshold * 100, color="red", linestyle="--", label=f"Threshold: {threshold * 100:g}%")
    ax.set(xlabel="BUSCO completeness (%)", ylabel="Number of samples")
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.legend()
    fig.tight_layout()
    fig.savefig(out / "busco_completeness.svg")
    plt.close(fig)
    report = {"created_at": now(), "input_runs": len(meta), "input_species": meta["species"].nunique(),
              "selected_runs": len(selected), "selected_species": selected["species"].nunique(),
              "busco_threshold": threshold, "requested_species": requested, "unknown_taxids": unknown,
              "metadata": file_record(metadata), "busco": file_record(busco),
              "taxonomy": file_record(taxonomy_db)}
    write_json(out / "selection.json", report)
    print(json.dumps({k: report[k] for k in ["input_runs", "selected_runs", "selected_species"]}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ["metadata", "busco", "cds-dir", "quant-dir", "taxonomy-db", "outdir"]:
        parser.add_argument(f"--{flag}", required=True)
    parser.add_argument("--threshold", type=float, default=0.5)
    parser.add_argument("--species-list")
    parser.add_argument("--missing-taxonomy", choices=["error", "allow"], default="error")
    prepare(**vars(parser.parse_args()))
