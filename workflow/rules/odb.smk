rule translate_cds:
    input:
        cds=lambda wc: species_row(wc)["cds"],
        code=f"{SCRIPTS}/translate_cds.py",
        common=f"{SCRIPTS}/common.py",
        cache_code=[f"{SCRIPTS}/protein_cache.py", f"{SCRIPTS}/dataset_assets.py"]
    output:
        protein=f"{PROTEINS}/{{species}}_protein.fa",
        provenance=f"{PROTEINS}/{{species}}_protein.json"
    params:
        seqkit="seqkit",
        table=config["translation"]["table"],
        cache=config["translation_cache"]
    threads: 1
    log: f"{LOG}/translate/{{species}}.log"
    conda: "../envs/seqkit.yaml"
    resources: mem_mb=2000
    shell:
        "{PYTHON:q} {input.code:q} --cds {input.cds:q} --output {output.protein:q} "
        "--provenance {output.provenance:q} --seqkit {params.seqkit:q} --table {params.table} "
        "--threads {threads} --cache-dir {params.cache:q} > {log:q} 2>&1"


rule make_manifests:
    input:
        samples=f"{META}/samples.tsv",
        code=f"{SCRIPTS}/make_manifests.py",
        common=f"{SCRIPTS}/common.py"
    output: directory(MANIFESTS)
    params: proteins=PROTEINS, chunk_size=config["odb"].get("chunk_size", 100)  # Species per chunk; also used by DAG planning.
    log: f"{LOG}/{ORTHOGROUP_MAPPING}/manifests.log"
    conda: "../envs/analysis.yaml"
    resources: mem_mb=1000
    shell:
        "{PYTHON:q} {input.code:q} --samples {input.samples:q} --protein-dir {params.proteins:q} "
        "--outdir {output:q} --chunk-size {params.chunk_size} > {log:q} 2>&1"


rule prepare_odb_reference:
    input:
        code=f"{SCRIPTS}/prepare_odb_reference.py",
        helpers=[f"{SCRIPTS}/common.py", f"{SCRIPTS}/odb_environment.py"]
    output:
        reference=f"{REFERENCE}/reference.json",
        inventory=f"{REFERENCE}/files.json",
        info=f"{REFERENCE}/dbinfo.txt",
        settings=f"{REFERENCE}/config.txt"
    params:
        root=REFERENCE,
        command="ODB-mapper",
        prefix="",
        version=ODB_VERSION,
        node=config["odb"]["node"]
    log: f"{LOG}/{ORTHOGROUP_MAPPING}/reference.log"
    conda: "../envs/odb.yaml"
    threads: 1
    resources: mem_mb=32000
    shell:
        "{PYTHON:q} {input.code:q} --reference-dir {params.root:q} --command {params.command:q} "
        "--prefix={params.prefix:q} --version {params.version:q} --node {params.node} > {log:q} 2>&1"


rule odb_map:
    input:
        manifests=odb_manifests,
        proteins=lambda wc: chunk_row(wc)["proteins"],
        reference=f"{REFERENCE}/reference.json",
        inventory=f"{REFERENCE}/files.json",
        code=f"{SCRIPTS}/run_odb_chunk.py",
        helpers=[f"{SCRIPTS}/odb_map.sh", f"{SCRIPTS}/odb_environment.py", f"{SCRIPTS}/common.py", f"{SCRIPTS}/make_manifests.py", f"{SCRIPTS}/incremental_odb.py"],
        samples=f"{META}/samples.tsv"
    output:
        annotations=f"{CHUNKS}/{{chunk}}/{{chunk}}.og.annotations",
        hits=f"{CHUNKS}/{{chunk}}/{{chunk}}.og.hits",
        summary=f"{CHUNKS}/{{chunk}}/{{chunk}}.summary.txt",
        provenance=f"{CHUNKS}/{{chunk}}/provenance.json"
    params:
        manifest=lambda wc: f"{ODB_PLAN if INCREMENTAL_ODB else MANIFESTS}/{wc.chunk}.fs",
        out=lambda wc: f"{CHUNKS}/{wc.chunk}",
        work=f"{WORK}/{ORTHOGROUP_MAPPING}",
        command="ODB-mapper",
        prefix="",
        version=ODB_VERSION,
        node=config["odb"]["node"],
        batch=lambda wildcards, threads: 4 * threads,  # Internal jobs per batch.
        incremental=int(INCREMENTAL_ODB),
        publish=([PYTHON, f"{SCRIPTS}/incremental_odb.py", "publish", "--samples", f"{META}/samples.tsv",
                  "--protein-dir", PROTEINS, "--cache-dir", ODB_CACHE, "--version", ODB_VERSION,
                  "--node", str(config["odb"]["node"])] if INCREMENTAL_ODB else [])
    threads: 16
    resources: mem_mb=192000
    log: f"{LOG}/{ORTHOGROUP_MAPPING}/chunks/{{chunk}}.log"
    benchmark: f"{LOG}/{ORTHOGROUP_MAPPING}/benchmarks/{{chunk}}.tsv"
    conda: "../envs/odb.yaml"
    shell:
        "{PYTHON:q} {input.code:q} --manifest {params.manifest:q} --reference {input.reference:q} "
        "--output-dir {params.out:q} --work-dir {params.work:q} --label {wildcards.chunk:q} "
        "--command {params.command:q} --prefix={params.prefix:q} --version {params.version:q} "
        "--node {params.node} --jobs {threads} --batch-size {params.batch} > {log:q} 2>&1; "
        "if [ {params.incremental} -eq 1 ]; then "
        "{params.publish:q} --chunk-dir {params.out:q} --label {wildcards.chunk:q} >> {log:q} 2>&1; fi"


rule collect_odb:
    input:
        samples=f"{META}/samples.tsv",
        manifests=odb_manifests,
        annotations=lambda wc: odb_source_files(wc, "og.annotations"),
        hits=lambda wc: odb_source_files(wc, "og.hits"),
        summaries=lambda wc: odb_source_files(wc, "summary.txt"),
        provenance=lambda wc: odb_source_files(wc, "provenance.json"),
        proteins=lambda wc: sorted({f'{PROTEINS}/{r["odb_species"]}_protein.fa' for r in sample_rows(wc)}),
        code=f"{SCRIPTS}/mapping_tables.py",
        helpers=[f"{SCRIPTS}/common.py", f"{SCRIPTS}/translate_cds.py", f"{SCRIPTS}/incremental_odb.py", f"{SCRIPTS}/dataset_assets.py"]
    output:
        snapshot=f"{MAPPING}/snapshot.json",
        tables=directory(f"{MAPPING}/species")
    params:
        out=MAPPING, cache=ODB_CACHE, chunks=CHUNKS, proteins=PROTEINS, plan=f"{MANIFESTS}/chunks.json",
        existing=(["--source-plan", f"{ODB_PLAN}/plan.json"] if INCREMENTAL_ODB else
                  ["--existing", EXISTING_ODB] if EXISTING_ODB else []),
        version=ODB_VERSION, node=config["odb"]["node"]
    log: f"{LOG}/{ORTHOGROUP_MAPPING}/tables.log"
    conda: "../envs/analysis.yaml"
    threads: 8
    resources: mem_mb=16000
    shell:
        "{PYTHON:q} {input.code:q} --samples {input.samples:q} --chunks {params.plan:q} "
        "--chunk-dir {params.chunks:q} --protein-dir {params.proteins:q} --outdir {params.out:q} --cache-dir {params.cache:q} --threads {threads} "
        "{params.existing:q} "
        "--version {params.version:q} --node {params.node} > {log:q} 2>&1"


checkpoint plan_incremental_odb:
    input:
        samples=f"{META}/samples.tsv",
        proteins=lambda wc: sorted({f'{PROTEINS}/{r["odb_species"]}_protein.fa' for r in sample_rows(wc)}),
        snapshots=[f"{EXISTING_ODB}/snapshot.json", f"{EXISTING_ODB}/annotations.tsv"] if EXISTING_ODB else [],
        code=f"{SCRIPTS}/incremental_odb.py",
        helpers=[f"{SCRIPTS}/common.py", f"{SCRIPTS}/make_manifests.py"]
    output: directory(ODB_PLAN)
    params:
        proteins=PROTEINS, cache=ODB_CACHE, reference=f"{REFERENCE}/reference.json",
        existing=["--existing", EXISTING_ODB] if EXISTING_ODB else [],
        version=ODB_VERSION, node=config["odb"]["node"], chunk_size=config["odb"].get("chunk_size", 20)
    log: f"{LOG}/{ORTHOGROUP_MAPPING}/incremental_plan.log"
    conda: "../envs/analysis.yaml"
    resources: mem_mb=2000
    shell:
        "{PYTHON:q} {input.code:q} plan --samples {input.samples:q} --protein-dir {params.proteins:q} "
        "--outdir {output:q} --cache-dir {params.cache:q} --reference {params.reference:q} "
        "--version {params.version:q} --node {params.node} --chunk-size {params.chunk_size} "
        "{params.existing:q} > {log:q} 2>&1"
