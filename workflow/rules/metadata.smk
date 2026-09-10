# Existing snapshots are external inputs, even if code or source settings change.
# Only register a producer when the fixed snapshot is missing.
if not Path(TAXONOMY_DB).exists():
    rule prepare_taxonomy:
        input:
            source=[config["taxonomy"]["source"]] if config["taxonomy"].get("source") else [],
            code=f"{SCRIPTS}/prepare_taxonomy.py",
            helpers=[f"{SCRIPTS}/snapshot_taxonomy.py", f"{SCRIPTS}/common.py"]
        output:
            database=TAXONOMY_DB,
            provenance=f"{TAXONOMY_DB}.json"
        params:
            source_flag="--source" if config["taxonomy"].get("source") else ""
        log: f"{LOG}/taxonomy_reference.log"
        conda: "../envs/analysis.yaml"
        threads: 1
        resources: mem_mb=8000
        shell:
            "{PYTHON:q} {input.code:q} --destination {output.database:q} "
            "{params.source_flag} {input.source:q} > {log:q} 2>&1"


checkpoint select_metadata:
    input:
        metadata=config["inputs"]["metadata"],
        busco=config["inputs"]["busco"],
        taxonomy=TAXONOMY_DB,
        subset=[config["selection"]["species_list"]] if config["selection"]["species_list"] else [],
        code=f"{SCRIPTS}/prepare_metadata.py",
        common=f"{SCRIPTS}/common.py"
    output:
        samples=f"{META}/samples.tsv",
        all_metadata=f"{META}/metadata_all.tsv",
        selected=f"{META}/metadata_high_busco.tsv",
        species=f"{META}/species_high_busco.txt",
        qc=f"{META}/selection.json",
        plot=f"{META}/busco_completeness.svg"
    params:
        outdir=META,
        cds=config["inputs"]["cds_dir"],
        quant=config["inputs"]["quant_dir"],
        threshold=config["selection"]["busco_threshold"],
        missing_taxonomy=config["selection"]["missing_taxonomy"],
        subset_flag="--species-list" if config["selection"]["species_list"] else ""
    log: f"{LOG}/metadata.log"
    conda: "../envs/analysis.yaml"
    resources: mem_mb=4000
    shell:
        "{PYTHON:q} {input.code:q} --metadata {input.metadata:q} --busco {input.busco:q} "
        "--taxonomy-db {input.taxonomy:q} --cds-dir {params.cds:q} --quant-dir {params.quant:q} "
        "--outdir {params.outdir:q} --threshold {params.threshold} --missing-taxonomy {params.missing_taxonomy:q} "
        "{params.subset_flag} {input.subset:q} > {log:q} 2>&1"


rule record_run:
    input:
        selection=f"{META}/selection.json",
        code=f"{SCRIPTS}/record_run.py",
        sources=sorted(str(p) for p in Path(workflow.basedir).rglob("*") if p.is_file()
                       and (p.suffix in {".py", ".sh", ".smk", ".yaml"} or p.name == "Snakefile"))
    output: f"{OUT}/run.json"
    params:
        resolved=json.dumps(config, sort_keys=True),
        workflow_dir=workflow.basedir
    log: f"{LOG}/provenance.log"
    conda: "../envs/analysis.yaml"
    resources: mem_mb=1000
    shell:
        "{PYTHON:q} {input.code:q} --config-json {params.resolved:q} --selection {input.selection:q} "
        "--workflow-dir {params.workflow_dir:q} --output {output:q} > {log:q} 2>&1"
