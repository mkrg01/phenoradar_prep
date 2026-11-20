"""Canonical names for trait clades, contrast pairs, and phylogenetic blocks."""

from pathlib import Path

from common import read_tsv

BLOCK_COLUMN = "phylogenetic_block_id"
LEGACY_BLOCK_COLUMN = "expanded_pair_clade_id"
CLADE_COLUMN = "trait_clade_id"


def canonical_column(name, *, trait_clades=False):
    if name == LEGACY_BLOCK_COLUMN or name.startswith(LEGACY_BLOCK_COLUMN + "_"):
        return BLOCK_COLUMN + name[len(LEGACY_BLOCK_COLUMN):]
    if trait_clades:
        return {"group": CLADE_COLUMN,
                "n_species_in_group": "n_species_in_trait_clade",
                "n_samples_in_group": "n_samples_in_trait_clade"}.get(name, name)
    return name


def canonical_row(row, *, trait_clades=False):
    """Read legacy names without modifying saved tables; reject conflicting aliases."""
    result = {}
    for name, value in row.items():
        canonical = canonical_column(name, trait_clades=trait_clades)
        if canonical in result:
            current = "" if result[canonical] is None else str(result[canonical]).strip()
            legacy = "" if value is None else str(value).strip()
            if current != legacy:
                raise ValueError(f"conflicting metadata aliases for {canonical}")
        result[canonical] = value
    return result


def canonical_rows(rows, *, trait_clades=False):
    return [canonical_row(row, trait_clades=trait_clades) for row in rows]


def read_blocks(directory):
    """Accept the former block-table filename in saved results."""
    directory = Path(directory)
    current = directory / "phylogenetic_blocks.tsv"
    legacy = directory / "expanded_pair_clades.tsv"
    rows = canonical_rows(read_tsv(current if current.is_file() else legacy))
    if current.is_file() and legacy.is_file() and rows != canonical_rows(read_tsv(legacy)):
        raise ValueError("conflicting phylogenetic block tables")
    return rows
