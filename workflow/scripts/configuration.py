"""Accepted configuration keys; defaults live in config/config.yaml."""

KEYS = {
    "": "analysis inputs selection taxonomy translation odb tpm alignment kegg phylogeny "
        "contrast taxonomy_audit exclude_species phenoradar",
    "inputs": "metadata species_trait busco cds_dir quant_dir",
    "selection": "busco_threshold species_list missing_taxonomy",
    "taxonomy": "source",
    "translation": "table",
    "odb": "existing_results node chunk_size threads batch_size mem_gb min_free_gb "
        "reference_min_free_gb allow_nonlocal keep_work",
    "tpm": "multimap",
    "alignment": "enabled threads mem_gb",
    "kegg": "enabled threads mem_gb ambiguity",
    "phylogeny": "enabled species_sets trait busco_full_dir busco_full_suffix lineage sequence_mode "
        "sequence_dir sequence_suffix outgroup max_markers min_protein_length max_unknown_fraction "
        "min_taxa trimal_mode align_threads tree_threads astral_threads preparation_mem_gb "
        "alignment_mem_gb trimming_mem_gb tree_mem_gb astral_mem_gb seed dating",
    "phylogeny.dating": "enabled calibration_source calibrations mem_gb lsd2 timetree",
    "phylogeny.dating.lsd2": "variance variance_parameter numsites",
    "phylogeny.dating.timetree": "representatives max_representatives max_queries min_studies offline",
    "contrast": "enabled trait",
    "taxonomy_audit": "enabled ranks outlierlevel collapse_monophyletic mem_gb",
    "phenoradar": "orthogroups kegg alignments kegg_groups orthogroup_annotations contrast tree",
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
