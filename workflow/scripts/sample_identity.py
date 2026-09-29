"""Separate biological species from independently assembled RNA-seq samples.

The historical ``species`` column in computational tables is a sample key.
Never recover taxonomy by splitting that key: accessions may contain underscores.
"""
import re

SAFE = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]*\Z")


def species_id(row):
    value = row.get("species_id") or row.get("scientific_name", "").replace(" ", "_")
    return value or row["species"]


def sample_id(name, run):
    taxon = name.replace(" ", "_")
    if not SAFE.fullmatch(taxon) or not SAFE.fullmatch(run):
        raise ValueError(f"unsafe species/run: {name!r}, {run!r}")
    return taxon + "_" + run


def annotate(row):
    """Validate explicit identities and return a copy with canonical columns."""
    row = dict(row)
    taxon = row["scientific_name"].replace(" ", "_")
    sample = sample_id(row["scientific_name"], row["run"])
    for key, value in (("species_id", taxon), ("analysis_sample_id", sample)):
        if row.get(key) and row[key] != value:
            raise ValueError(f"{key} differs from scientific_name/run: {row[key]}")
        row[key] = value
    return row


def traits_for_samples(rows, annotations):
    return {r["species"]: annotations[species_id(r)] for r in rows
            if species_id(r) in annotations}


def select_samples(rows, requested):
    """Species selectors include all their samples; sample selectors are exact."""
    rows = list(rows)
    aliases = {}
    for row in rows:
        name = row["species"]
        for key in {name, species_id(row)}:
            aliases.setdefault(key, set()).add(name)
    unknown = set(requested) - aliases.keys()
    if unknown:
        raise ValueError(f"unknown species/sample IDs: {sorted(unknown)}")
    return set().union(*(aliases[n] for n in requested)) if requested else set()
