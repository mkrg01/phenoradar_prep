rule translate_cds:
    input:
        cds=lambda wc: species_row(wc)["cds"],
        code=f"{SCRIPTS}/translate_cds.py",
        common=f"{SCRIPTS}/common.py"
    output:
        protein=f"{PROTEINS}/{{species}}_protein.fa",
        provenance=f"{PROTEINS}/{{species}}_protein.json"
    params:
        seqkit=config["tools"]["seqkit"],
        table=config["translation"]["table"]
    threads: 1
    log: f"{LOG}/translate/{{species}}.log"
    conda: "../envs/seqkit.yaml"
    resources: mem_mb=2000
    shell:
        "{PYTHON:q} {input.code:q} --cds {input.cds:q} --output {output.protein:q} "
        "--provenance {output.provenance:q} --seqkit {params.seqkit:q} --table {params.table} "
        "--threads {threads} > {log:q} 2>&1"


rule make_manifests:
    input:
        samples=f"{META}/samples.tsv",
        code=f"{SCRIPTS}/make_manifests.py",
        common=f"{SCRIPTS}/common.py"
    output: directory(MANIFESTS)
    params: proteins=PROTEINS, chunk_size=config["odb"]["chunk_size"]
    log: f"{LOG}/manifests.log"
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
        command=config["tools"]["odb_command"],
        prefix=config["tools"]["odb_prefix"],
        version=config["odb"]["version"],
        node=config["odb"]["node"],
        free=config["odb"]["reference_min_free_gb"],
        storage="--allow-nonlocal" if config["odb"]["allow_nonlocal"] else ""
    log: f"{LOG}/odb_reference.log"
    conda: "../envs/odb.yaml"
    threads: 1
    resources: mem_mb=32000, odb_slots=1
    shell:
        "{PYTHON:q} {input.code:q} --reference-dir {params.root:q} --command {params.command:q} "
        "--prefix={params.prefix:q} --version {params.version:q} --node {params.node} "
        "--min-free-gb {params.free} {params.storage} > {log:q} 2>&1"


rule odb_map:
    input:
        manifests=MANIFESTS,
        proteins=lambda wc: chunk_row(wc)["proteins"],
        reference=f"{REFERENCE}/reference.json",
        inventory=f"{REFERENCE}/files.json",
        code=f"{SCRIPTS}/run_odb_chunk.py",
        helpers=[f"{SCRIPTS}/odb_map.sh", f"{SCRIPTS}/odb_environment.py", f"{SCRIPTS}/common.py", f"{SCRIPTS}/make_manifests.py"]
    output:
        annotations=f"{CHUNKS}/{{chunk}}/{{chunk}}.og.annotations",
        hits=f"{CHUNKS}/{{chunk}}/{{chunk}}.og.hits",
        summary=f"{CHUNKS}/{{chunk}}/{{chunk}}.summary.txt",
        provenance=f"{CHUNKS}/{{chunk}}/provenance.json"
    params:
        manifest=lambda wc: f"{MANIFESTS}/{wc.chunk}.fs",
        out=lambda wc: f"{CHUNKS}/{wc.chunk}",
        work=f"{WORK}/odb",
        command=config["tools"]["odb_command"],
        prefix=config["tools"]["odb_prefix"],
        version=config["odb"]["version"],
        node=config["odb"]["node"],
        batch=config["odb"]["batch_size"],
        free=config["odb"]["min_free_gb"],
        storage="--allow-nonlocal" if config["odb"]["allow_nonlocal"] else "",
        keep="--keep-work" if config["odb"]["keep_work"] else ""
    threads: config["odb"]["threads"]
    resources:
        mem_mb=config["odb"]["mem_gb"] * 1000,
        odb_slots=1
    log: f"{LOG}/odb/{{chunk}}.log"
    benchmark: f"{LOG}/benchmarks/odb_{{chunk}}.tsv"
    conda: "../envs/odb.yaml"
    shell:
        "{PYTHON:q} {input.code:q} --manifest {params.manifest:q} --reference {input.reference:q} "
        "--output-dir {params.out:q} --work-dir {params.work:q} --label {wildcards.chunk:q} "
        "--command {params.command:q} --prefix={params.prefix:q} --version {params.version:q} "
        "--node {params.node} --jobs {threads} --batch-size {params.batch} --min-free-gb {params.free} "
        "{params.storage} {params.keep} > {log:q} 2>&1"


rule merge_odb:
    input:
        samples=f"{META}/samples.tsv",
        manifests=MANIFESTS,
        annotations=lambda wc: [f"{EXISTING_ODB}/annotations.tsv"] if EXISTING_ODB else chunk_outputs(wc, "og.annotations"),
        hits=lambda wc: [] if EXISTING_ODB else chunk_outputs(wc, "og.hits"),
        summaries=lambda wc: [] if EXISTING_ODB else chunk_outputs(wc, "summary.txt"),
        provenance=lambda wc: [f"{EXISTING_ODB}/snapshot.json"] if EXISTING_ODB else
            [f'{CHUNKS}/{r["chunk"]}/provenance.json' for r in chunk_rows(wc)],
        proteins=lambda wc: sorted({f'{PROTEINS}/{r["odb_species"]}_protein.fa' for r in sample_rows(wc)}),
        code=f"{SCRIPTS}/merge_odb.py",
        helpers=[f"{SCRIPTS}/common.py", f"{SCRIPTS}/translate_cds.py"]
    output:
        database=f"{MERGED}/mappings.sqlite",
        mappings=f"{MERGED}/gene_orthogroups.tsv",
        qc=f"{MERGED}/merge_qc.json"
    params:
        chunks=CHUNKS, proteins=PROTEINS, plan=f"{MANIFESTS}/chunks.json",
        existing=["--existing", EXISTING_ODB] if EXISTING_ODB else [],
        version=config["odb"]["version"], node=config["odb"]["node"]
    log: f"{LOG}/merge_odb.log"
    conda: "../envs/analysis.yaml"
    resources: mem_mb=8000
    shell:
        "{PYTHON:q} {input.code:q} --samples {input.samples:q} --chunks {params.plan:q} "
        "--chunk-dir {params.chunks:q} --protein-dir {params.proteins:q} --database {output.database:q} "
        "--mappings {output.mappings:q} --qc {output.qc:q} {params.existing:q} "
        "--version {params.version:q} --node {params.node} > {log:q} 2>&1"
