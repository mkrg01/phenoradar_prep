# Existing snapshots are shared external inputs, even if code or run_name changes.
# Only register a producer when the fixed snapshot is missing.
if not Path(TAXONOMY_DB).exists():
    rule prepare_taxonomy:
        input:
            code=f"{SCRIPTS}/prepare_taxonomy.py",
            helpers=[f"{SCRIPTS}/snapshot_taxonomy.py", f"{SCRIPTS}/common.py", f"{SCRIPTS}/sample_identity.py"]
        output:
            database=TAXONOMY_DB,
            provenance=f"{TAXONOMY_DB}.json"
        log: f"{LOG}/taxonomy_reference.log"
        conda: "../envs/analysis.yaml"
        threads: 1
        resources: mem_mb=8000
        shell:
            "{PYTHON:q} {input.code:q} --destination {output.database:q} "
            "> {log:q} 2>&1"


checkpoint select_metadata:
    input:
        metadata=INPUTS["metadata"],
        busco=INPUTS["busco"],
        taxonomy=TAXONOMY_DB,
        subset=[SPECIES_LIST] if config["selection"]["species_list"] else [],
        code=f"{SCRIPTS}/prepare_metadata.py",
        common=[f"{SCRIPTS}/common.py", f"{SCRIPTS}/sample_identity.py"]
    output:
        samples=f"{META}/samples.tsv",
        all_metadata=f"{META}/metadata_all.tsv",
        selected=f"{META}/metadata_high_busco.tsv",
        species=f"{META}/species_high_busco.txt",
        qc=f"{META}/selection.json",
        plot=f"{META}/busco_completeness.svg"
    params:
        outdir=META,
        cds=INPUTS["cds_dir"],
        quant=INPUTS["quant_dir"],
        threshold=config["selection"]["busco_threshold"],
        subset_flag="--species-list" if config["selection"]["species_list"] else ""
    log: f"{LOG}/metadata.log"
    conda: "../envs/analysis.yaml"
    resources: mem_mb=4000
    shell:
        "{PYTHON:q} {input.code:q} --metadata {input.metadata:q} --busco {input.busco:q} "
        "--taxonomy-db {input.taxonomy:q} --cds-dir {params.cds:q} --quant-dir {params.quant:q} "
        "--outdir {params.outdir:q} --threshold {params.threshold} "
        "{params.subset_flag} {input.subset:q} > {log:q} 2>&1"


# Missing optional traits leave base metadata blank; adding the file creates a dependency.
PHENORADAR_TRAITS = INPUTS["species_trait"] if Path(INPUTS["species_trait"]).is_file() else ""


rule normalize_traits:
    input:
        samples=f"{META}/samples.tsv",
        traits=[PHENORADAR_TRAITS] if PHENORADAR_TRAITS else [],
        code=f"{SCRIPTS}/trait_encoding.py",
        helpers=[f"{SCRIPTS}/{name}.py" for name in ["species_traits", "common", "sample_identity"]]
    output:
        traits=f"{META}/species_trait.tsv",
        encoding=f"{META}/trait_encoding.json"
    params:
        trait=config["trait"], mapping=json.dumps(config.get("trait_encoding")),
        traits_flag="--traits" if PHENORADAR_TRAITS else ""
    log: f"{LOG}/trait_encoding.log"
    conda: "../envs/analysis.yaml"
    resources: mem_mb=1000
    shell:
        "{PYTHON:q} {input.code:q} --samples {input.samples:q} --output {output.traits:q} "
        "--report {output.encoding:q} --trait {params.trait:q} --mapping {params.mapping:q} "
        "{params.traits_flag} {input.traits:q} > {log:q} 2>&1"


rule prepare_phenoradar_metadata:
    input:
        samples=f"{META}/samples.tsv",
        metadata=f"{META}/metadata_high_busco.tsv",
        traits=rules.normalize_traits.output.traits,
        encoding=rules.normalize_traits.output.encoding,
        code=f"{SCRIPTS}/phenoradar_metadata.py",
        helpers=[f"{SCRIPTS}/species_traits.py", f"{SCRIPTS}/common.py", f"{SCRIPTS}/sample_identity.py"]
    output: f"{OUT}/{BASE_METADATA}"
    params: trait=config["trait"]
    log: f"{LOG}/phenoradar_metadata.log"
    conda: "../envs/analysis.yaml"
    resources: mem_mb=1000
    shell:
        "{PYTHON:q} {input.code:q} --samples {input.samples:q} --metadata {input.metadata:q} "
        "--output {output:q} --trait {params.trait:q} "
        "--traits {input.traits:q} > {log:q} 2>&1"


def consumer_pair_tables():
    # Only configured branches are inputs; stale contrast files are ignored.
    return ({branch: f"{OUT}/{TREE_BRANCHES[branch]}/contrast/species_metadata.tsv"
             for branch in PHY["trees"]} if PHY["contrast_pairs"]["enabled"] else {})


rule prepare_consumer_metadata:
    input:
        base=f"{OUT}/{BASE_METADATA}",
        taxonomy=f"{META}/metadata_high_busco.tsv",
        pairs=lambda wc: list(consumer_pair_tables().values()),
        code=f"{SCRIPTS}/consumer_metadata.py",
        helpers=[f"{SCRIPTS}/{name}.py" for name in ["common", "phenoradar_metadata", "sample_identity", "species_traits"]]
    output:
        table=f"{OUT}/{CONSUMER_METADATA}",
        provenance=f"{OUT}/{Path(CONSUMER_METADATA).with_suffix('.json')}"
    params: trait=config["trait"], pairs=lambda wc: json.dumps(consumer_pair_tables(), sort_keys=True)
    log: f"{LOG}/consumer_metadata.log"
    conda: "../envs/analysis.yaml"
    resources: mem_mb=1000
    shell:
        "{PYTHON:q} {input.code:q} --base {input.base:q} --taxonomy {input.taxonomy:q} "
        "--output {output.table:q} --trait {params.trait:q} --pair-tables {params.pairs:q} > {log:q} 2>&1"


rule record_run:
    input:
        selection=f"{META}/selection.json",
        code=f"{SCRIPTS}/record_run.py",
        sources=sorted(str(p) for p in Path(workflow.basedir).rglob("*") if p.is_file()
                       and (p.suffix in {".py", ".sh", ".smk", ".yaml"} or p.name == "Snakefile"))
    output: f"{OUT}/run.json"
    params:
        resolved=json.dumps(config, sort_keys=True),
        image=["--container-image", CONTAINER_IMAGE] if CONTAINER_IMAGE else [],
        workflow_dir=workflow.basedir
    log: f"{LOG}/provenance.log"
    conda: "../envs/analysis.yaml"
    resources: mem_mb=1000
    shell:
        "{PYTHON:q} {input.code:q} --config-json {params.resolved:q} --selection {input.selection:q} "
        "--workflow-dir {params.workflow_dir:q} {params.image:q} --output {output:q} > {log:q} 2>&1"
