"""The species-level TSV is the sole source of phenotype annotations."""
import argparse
import csv
from pathlib import Path
import re

from common import file_record, now, read_tsv, write_json, write_tsv


def read_species_traits(path, column="C4"):
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {"species", column} <= set(reader.fieldnames or []):
            raise ValueError(f"species trait table requires species and {column} columns")
        result = {}
        for row in reader:
            name = row["species"].strip().replace(" ", "_")
            if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", name):
                raise ValueError(f"invalid species trait label: {name}")
            if name in result:
                raise ValueError(f"duplicate species after normalization in trait table: {name}")
            value = row[column].strip()
            value = "" if value in {"", "NA", "NaN", "nan"} else value
            if column == "C4" and value not in {"", "0", "1"}:
                raise ValueError(f"C4 must be 0, 1 or missing: {name}={value}")
            result[name] = value
    if not result:
        raise ValueError("empty species trait table")
    return result


def select_phenotyped(samples, traits, outdir, trait="C4", min_taxa=4):
    """Keep every selected species with an observed trait, including state zero."""
    rows = read_tsv(samples)
    annotation = read_species_traits(traits, trait)
    species = {r["species"] for r in rows}
    eligible = {n for n in species if annotation.get(n, "") != ""}
    if len(eligible) < min_taxa:
        raise ValueError(f"phylogeny requires at least {min_taxa} phenotyped species; found {len(eligible)} for {trait}")
    out = Path(outdir)
    write_tsv(out / "samples.tsv", list(rows[0]), [r for r in rows if r["species"] in eligible])
    write_json(out / "selection.json", {
        "created_at": now(), "samples": file_record(samples), "species_trait": file_record(traits),
        "trait": trait, "min_taxa": min_taxa, "dataset_species": len(species),
        "inference_species": len(eligible), "missing_trait_species": sorted(species - eligible),
        "missing_trait_rows": sorted(species - set(annotation)),
    })


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["samples", "traits", "outdir"]:
        parser.add_argument(f"--{name}", required=True)
    parser.add_argument("--trait", default="C4")
    parser.add_argument("--min-taxa", type=int, default=4)
    select_phenotyped(**vars(parser.parse_args()))
