def phenoradar_snapshot():
    # Consume a completed snapshot. Absolute paths deliberately avoid scheduling
    # upstream producers, including optional KEGG, alignment and tree workflows.
    return discover_phenoradar_inputs(OUT, PHENORADAR, EXCLUDE_SPECIES)


rule phenoradar_inputs:
    input:
        snapshot=lambda wc: phenoradar_snapshot()["files"],
        code=f"{SCRIPTS}/phenoradar_inputs.py",
        helpers=[f"{SCRIPTS}/{name}.py" for name in
                 ["phenoradar_metadata", "species_traits", "filter_species", "layout", "common"]]
    # This explicit manual target always validates and refreshes its publication.
    # Do not declare directory() output: Snakemake would delete the previous
    # bundle before validation, and again on failure, bypassing atomic rollback.
    params:
        bundle=f"{OUT}/phenoradar_inputs",
        source=str(Path(OUT).resolve()),
        settings=json.dumps(PHENORADAR, sort_keys=True),
        excluded=json.dumps(EXCLUDE_SPECIES),
        trait=config["contrast"]["trait"],
        inventory=lambda wc: json.dumps(phenoradar_snapshot()["sections"], sort_keys=True)
    conda: "../envs/analysis.yaml"
    resources: mem_mb=4000
    log: f"{LOG}/phenoradar_inputs.log"
    shell:
        "{PYTHON:q} {input.code:q} --source {params.source:q} --outdir {params.bundle:q} "
        "--settings {params.settings:q} --exclude-species {params.excluded:q} "
        "--trait {params.trait:q} > {log:q} 2>&1"
