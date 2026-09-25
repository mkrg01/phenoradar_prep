"""Accepted configuration keys, including optional tool-specific overrides."""
KEYS = {
    "": "run_name selection translation odb tpm alignment kegg phylogeny "
        "trait exclude_species seed input_root build_manifest translation_cache",
    "selection": "busco_threshold species_list missing_taxonomy",
    "translation": "table",
    "odb": "existing_results node incremental cache_dir chunk_size",
    "tpm": "multimap",
    "alignment": "enabled",
    "kegg": "enabled ambiguity",
    "phylogeny": "trees lineage "
        "outgroup max_markers min_protein_length max_unknown_fraction "
        "min_taxa trimal_mode dating contrast_pairs taxonomy_check",
    "phylogeny.dating": "enabled calibration_source treepl",
    "phylogeny.dating.treepl": "smooth cvstart cvstop cvmultstep lfiter pliter cviter thorough",
    "phylogeny.contrast_pairs": "enabled",
    "phylogeny.taxonomy_check": "enabled ranks outlierlevel collapse_monophyletic",
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
    existing = config.get("odb", {}).get("existing_results")
    if existing is not None and (not isinstance(existing, str) or not existing.strip()):
        raise ValueError("odb.existing_results must be null or a snapshot directory")
    if "input_root" in config and (not isinstance(config["input_root"], str) or not config["input_root"].strip()):
        raise ValueError("input_root must be a directory")
    odb = config.get("odb", {})
    if type(odb.get("incremental", False)) is not bool:
        raise ValueError("odb.incremental must be true or false")
    if type(odb.get("chunk_size", 100)) is not int or odb.get("chunk_size", 100) < 1:
        raise ValueError("odb.chunk_size must be a positive integer")
    if not isinstance(odb.get("cache_dir", "resources/odb_cache"), str) or not odb.get("cache_dir", "resources/odb_cache").strip():
        raise ValueError("odb.cache_dir must be a directory")
    # Use the positive signed 32-bit range supported by all seeded tools.
    if "seed" in config and (type(config["seed"]) is not int or not 1 <= config["seed"] <= 2147483647):
        raise ValueError("seed must be an integer between 1 and 2147483647")
    if "trait" in config:
        trait = config["trait"]
        if not isinstance(trait, str) or not trait.strip():
            raise ValueError("trait must name a column in input/species_trait.tsv")
        if trait in {"species", "role", "group", "representative", "is_representative",
                     "n_species_in_group", "contrast_pair_id"}:
            raise ValueError("trait conflicts with a reserved output column")
    phylogeny = config.get("phylogeny", {})
    if "trees" in phylogeny:
        trees = phylogeny["trees"]
        if (not isinstance(trees, list)
                or any(not isinstance(t, str) or t not in {"all", "phenotyped", "representatives"} for t in trees)
                or len(trees) != len(set(trees))):
            raise ValueError("phylogeny.trees must be a list of all, phenotyped, representatives, without duplicates; [] disables inference")
    for name in ["contrast_pairs", "dating", "taxonomy_check"]:
        if "enabled" in phylogeny.get(name, {}) and type(phylogeny[name]["enabled"]) is not bool:
            raise ValueError(f"phylogeny.{name}.enabled must be true or false")


def validate_analysis(config, targets=()):
    """Validate the resolved analysis and explicitly requested public targets."""
    phylogeny = config["phylogeny"]
    trees = phylogeny["trees"]
    for name in ["contrast_pairs", "dating", "taxonomy_check"]:
        if phylogeny[name]["enabled"] and not trees:
            raise ValueError(f"phylogeny.{name}.enabled requires a nonempty phylogeny.trees")
    if "representatives" in trees:
        for name in ["dating", "taxonomy_check"]:
            if phylogeny[name]["enabled"]:
                raise ValueError(f"phylogeny.{name} does not yet support representatives; use a separate all/phenotyped analysis")
    for target in targets:
        if target in {"phylogeny", "phylogeny_prepare", "phylogeny_calibrations", "timetree", "taxonomy_check", "contrast_pairs"}:
            if not trees:
                raise ValueError(f"target {target} requires a nonempty phylogeny.trees")
            section = {"contrast_pairs": "contrast_pairs", "timetree": "dating",
                       "phylogeny_calibrations": "dating", "taxonomy_check": "taxonomy_check"}.get(target)
            if section and not phylogeny[section]["enabled"]:
                raise ValueError(f"target {target} requires phylogeny.{section}.enabled: true; targets do not override configuration")
