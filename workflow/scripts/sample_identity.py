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


def sample_run_map(rows):
    """Resolve runs from explicit metadata, requiring one run per sample key."""
    samples, runs = {}, set()
    for row in rows:
        sample, run = row["species"], row["run"]
        if not sample or not run or run in runs:
            raise ValueError("empty sample/run or duplicate run in sample manifest")
        if sample in samples:
            raise ValueError("multiple runs per sample key; assign a distinct sample ID to each run: " + sample)
        samples[sample] = run
        runs.add(run)
    return samples


def expression_run(row, samples):
    """Resolve native rows and cross-check optional legacy run columns."""
    sample = row["species"]
    if sample not in samples:
        raise ValueError("expression sample is absent from sample manifest: " + sample)
    run = samples[sample]
    if "run" in row and row["run"] != run:
        raise ValueError("expression run/species differs from sample manifest: " + sample)
    return run


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
