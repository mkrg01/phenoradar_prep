def filter_snapshot():
    # Absolute input paths intentionally describe an already completed snapshot.
    # They do not request producers under results/<analysis>/, so this explicit
    # export cannot trigger old selection, mapping or inference checkpoints.
    return discover_filter_inputs(OUT, config["inputs"].get("species_trait"))


rule filter_species:
    input:
        snapshot=lambda wc: filter_snapshot()["files"],
        code=f"{SCRIPTS}/filter_species.py",
        phylogeny_code=f"{SCRIPTS}/filter_species_phylogeny.py",
        contrast_code=[f"{SCRIPTS}/{name}.py" for name in
                       ["contrast_pairs", "plot_contrast_tree", "species_traits", "phylogeny_root"]],
        common=f"{SCRIPTS}/common.py"
    output: bundle=directory(f"{OUT}/filtered")
    params:
        source=str(Path(OUT).resolve()),
        excluded=json.dumps(EXCLUDE_SPECIES),
        traits=config["inputs"].get("species_trait") or "",
        trait=config["contrast"]["trait"], seed=PHY["seed"],
        inventory=lambda wc: json.dumps(filter_snapshot()["sections"], sort_keys=True)
    conda: "../envs/timetree.yaml"
    resources: mem_mb=8000
    log: f"{LOG}/filter_species.log"
    shell:
        "{PYTHON:q} {input.code:q} --source {params.source:q} --outdir {output.bundle:q} "
        "--exclude-species {params.excluded:q} --traits {params.traits:q} "
        "--contrast-trait {params.trait:q} --seed {params.seed} > {log:q} 2>&1"
