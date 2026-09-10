"""The species-level TSV is the sole source of phenotype annotations."""
import csv
import re


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
