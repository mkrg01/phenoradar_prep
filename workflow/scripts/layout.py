"""Canonical input and result paths shared by workflow stages."""

INPUTS = {
    "metadata": "input/metadata.tsv",
    "species_trait": "input/species_trait.tsv",
    "busco": "input/busco/summary.tsv",
    "cds_dir": "input/cds",
    "quant_dir": "input/quant",
}
BUSCO_FULL = "input/busco/full"
SPECIES_LIST = "input/species_list.txt"
CALIBRATIONS = "input/calibrations.tsv"


def input_layout(root="input"):
    """Resolve the same input contract inside an immutable dataset snapshot."""
    from pathlib import Path
    root = Path(root)
    inputs = {key: str(root / Path(value).relative_to("input")) for key, value in INPUTS.items()}
    return inputs, str(root / "busco/full"), str(root / "species_list.txt"), str(root / "calibrations.tsv")

ORTHOGROUP_MAPPING = "orthogroups/mapping"
ORTHOGROUP_EXPRESSION = "orthogroups/expression"
ORTHOGROUP_ALIGNMENTS = "orthogroups/alignments"
PHYLOGENY_BRANCHES = {"all": "phylogeny/all", "phenotyped": "phylogeny/phenotyped"}
REPRESENTATIVES = "phylogeny/representatives"
ROOTING = f"{PHYLOGENY_BRANCHES['all']}/rooting"
CONTRAST_BRANCHES = frozenset(f"{branch}/contrast" for branch in
                             [*PHYLOGENY_BRANCHES.values(), REPRESENTATIVES])
