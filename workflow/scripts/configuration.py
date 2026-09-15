"""Accepted configuration keys, including optional tool-specific overrides."""
KEYS = {
    "": "run_name selection translation odb tpm alignment kegg phylogeny "
        "contrast taxonomy_check exclude_species seed",
    "selection": "busco_threshold species_list missing_taxonomy",
    "translation": "table",
    "odb": "node",
    "tpm": "multimap",
    "alignment": "enabled",
    "kegg": "enabled ambiguity",
    "phylogeny": "enabled species_sets trait lineage "
        "outgroup max_markers min_protein_length max_unknown_fraction "
        "min_taxa trimal_mode dating",
    "phylogeny.dating": "enabled calibration_source treepl",
    "phylogeny.dating.treepl": "smooth cvstart cvstop cvmultstep lfiter pliter cviter thorough",
    "contrast": "enabled trait",
    "taxonomy_check": "enabled ranks outlierlevel collapse_monophyletic",
}


def validate_keys(config):
    """Reject unknown keys and malformed sections before building the DAG."""
    for section, keys in KEYS.items():
        values = config
        for part in section.split(".") if section else []:
            values = values.get(part, {})
        if not isinstance(values, dict):
            raise ValueError(f"{section or 'configuration'} must be a mapping")
        unknown = sorted(set(values) - set(keys.split()), key=str)
        if unknown:
            paths = [f"{section}.{key}" if section else str(key) for key in unknown]
            raise ValueError("unknown configuration settings: " + ", ".join(paths))
    if "species_list" in config.get("selection", {}) and type(config["selection"]["species_list"]) is not bool:
        raise ValueError("selection.species_list must be true or false; true reads input/species_list.txt")
    # Use the positive signed 32-bit range supported by all seeded tools.
    if "seed" in config and (type(config["seed"]) is not int or not 1 <= config["seed"] <= 2147483647):
        raise ValueError("seed must be an integer between 1 and 2147483647")
