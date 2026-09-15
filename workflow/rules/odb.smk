rule translate_cds:
    input:
        cds=lambda wc: species_row(wc)["cds"],
        code=f"{SCRIPTS}/translate_cds.py",
        common=f"{SCRIPTS}/common.py"
    output:
        protein=f"{PROTEINS}/{{species}}_protein.fa",
        provenance=f"{PROTEINS}/{{species}}_protein.json"
    params:
        seqkit="seqkit",
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
    params: proteins=PROTEINS, chunk_size=100  # Species per chunk; also used by DAG planning.
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
        work=f"{WORK}/{ORTHOGROUP_MAPPING}",
        command="ODB-mapper",
        prefix="",
        version=ODB_VERSION,
        node=config["odb"]["node"],
        batch=lambda wildcards, threads: 4 * threads  # Internal jobs per batch.
    threads: 16
    resources: mem_mb=192000
    log: f"{LOG}/{ORTHOGROUP_MAPPING}/chunks/{{chunk}}.log"
    benchmark: f"{LOG}/{ORTHOGROUP_MAPPING}/benchmarks/{{chunk}}.tsv"
    conda: "../envs/odb.yaml"
    shell:
        "{PYTHON:q} {input.code:q} --manifest {params.manifest:q} --reference {input.reference:q} "
        "--output-dir {params.out:q} --work-dir {params.work:q} --label {wildcards.chunk:q} "
        "--command {params.command:q} --prefix={params.prefix:q} --version {params.version:q} "
        "--node {params.node} --jobs {threads} --batch-size {params.batch} > {log:q} 2>&1"


rule merge_odb:
    input:
        samples=f"{META}/samples.tsv",
        manifests=MANIFESTS,
        annotations=lambda wc: chunk_outputs(wc, "og.annotations"),
        hits=lambda wc: chunk_outputs(wc, "og.hits"),
        summaries=lambda wc: chunk_outputs(wc, "summary.txt"),
        provenance=lambda wc: [f'{CHUNKS}/{r["chunk"]}/provenance.json' for r in chunk_rows(wc)],
        proteins=lambda wc: sorted({f'{PROTEINS}/{r["odb_species"]}_protein.fa' for r in sample_rows(wc)}),
        code=f"{SCRIPTS}/merge_odb.py",
        helpers=[f"{SCRIPTS}/common.py", f"{SCRIPTS}/translate_cds.py"]
    output:
        database=f"{MAPPING}/mappings.sqlite",
        mappings=f"{MAPPING}/gene_orthogroups.tsv",
        qc=f"{MAPPING}/merge_qc.json"
    params:
        chunks=CHUNKS, proteins=PROTEINS, plan=f"{MANIFESTS}/chunks.json"
    log: f"{LOG}/{ORTHOGROUP_MAPPING}/merge.log"
    conda: "../envs/analysis.yaml"
    resources: mem_mb=8000
    shell:
        "{PYTHON:q} {input.code:q} --samples {input.samples:q} --chunks {params.plan:q} "
        "--chunk-dir {params.chunks:q} --protein-dir {params.proteins:q} --database {output.database:q} "
        "--mappings {output.mappings:q} --qc {output.qc:q} > {log:q} 2>&1"
