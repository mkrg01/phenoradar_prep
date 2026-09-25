"""Read manually curated run exclusions, also usable by metadata preparation tools."""
import csv
from pathlib import Path

from dataset_assets import SAFE


def read_exclusions(path):
    """Return exact run-accession decisions; accessions absent from metadata are valid."""
    if path is None:
        return {}
    with Path(path).open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames
        if not fields or "accession" not in fields or len(set(fields)) != len(fields) or any(not f.strip() for f in fields):
            raise ValueError("exclusion TSV requires unique columns including accession")
        result = {}
        for number, row in enumerate(reader, 2):
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"malformed exclusion TSV row: {path}:{number}")
            accession = row["accession"]
            if not SAFE.fullmatch(accession):
                raise ValueError(f"invalid excluded accession: {path}:{number}: {accession!r}")
            if accession in result:
                raise ValueError(f"duplicate excluded accession: {accession}")
            result[accession] = row
    return result


def partition(items, exclusions):
    """Filter metadata run IDs without guessing species/project-level exclusions."""
    included, excluded = [], []
    for item in items:
        run = item["row"]["run"]
        decision = exclusions.get(run)
        if decision is None:
            included.append(item)
        else:
            excluded.append({"species": item["species"], "run": run,
                             "reason": decision.get("reason", "") or "manually excluded"})
    return included, excluded
