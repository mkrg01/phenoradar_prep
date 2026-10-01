# A completed snapshot is an immutable input, independent of code changes.
# Failed downloads remain outside the snapshot and can be reused on retry.
if not Path(KEGG_REFERENCE).exists():
    rule prepare_kegg_reference:
        input:
            code=f"{SCRIPTS}/bootstrap_kegg_reference.py",
            helpers=[f"{SCRIPTS}/prepare_kegg_reference.py", f"{SCRIPTS}/verify_kegg_reference.py",
                     f"{SCRIPTS}/common.py"]
        output:
            reference=f"{KEGG_REFERENCE}/reference.json",
            inventory=f"{KEGG_REFERENCE}/files.json",
            modules=f"{KEGG_REFERENCE}/ko_modules.tsv",
            pathways=f"{KEGG_REFERENCE}/ko_pathways.tsv"
        params: root=KEGG_REFERENCE
        threads: 1
        resources: mem_mb=4000
        conda: "../envs/analysis.yaml"
        log: f"{LOG}/kegg/reference_prepare.log"
        shell:
            "{PYTHON:q} {input.code:q} --reference-dir {params.root:q} > {log:q} 2>&1"


rule kegg_verify_reference:
    input:
        reference=f"{KEGG_REFERENCE}/reference.json",
        inventory=f"{KEGG_REFERENCE}/files.json",
        modules=f"{KEGG_REFERENCE}/ko_modules.tsv",
        pathways=f"{KEGG_REFERENCE}/ko_pathways.tsv",
        code=f"{SCRIPTS}/publish_kegg_reference.py",
        verifier=f"{SCRIPTS}/verify_kegg_reference.py",
        preparation=f"{SCRIPTS}/prepare_kegg_reference.py",
        common=f"{SCRIPTS}/common.py"
    output:
        qc=f"{KEGG}/reference_qc.json",
        modules=f"{KEGG}/ko_modules.tsv",
        pathways=f"{KEGG}/ko_pathways.tsv"
    log: f"{LOG}/kegg/reference.log"
    conda: "../envs/analysis.yaml"
    resources: mem_mb=2000
    shell:
        "{PYTHON:q} {input.code:q} --reference {input.reference:q} --qc {output.qc:q} "
        "--modules {output.modules:q} --pathways {output.pathways:q} > {log:q} 2>&1"


checkpoint select_kegg_representatives:
    input:
        samples=f"{META}/samples.tsv",
        mapping=f"{MAPPING}/snapshot.json",
        tables=f"{MAPPING}/species",
        proteins=lambda wc: sorted({f'{PROTEINS}/{r["odb_species"]}_protein.fa' for r in sample_rows(wc)}),
        code=f"{SCRIPTS}/select_kegg_representatives.py",
        helpers=[f"{SCRIPTS}/{name}.py" for name in [
            "common", "mapping_tables", "dataset_assets", "busco_phylogeny", "aggregate_ko_tpm", "run_kofam"]]
    output: root=directory(KEGG_REPRESENTATIVES)
    params: proteins=PROTEINS, batch_size=5000
    resources: mem_mb=8000
    conda: "../envs/analysis.yaml"
    log: f"{LOG}/kegg/representatives.log"
    benchmark: f"{LOG}/kegg/benchmarks/representatives.tsv"
    shell:
        "{PYTHON:q} {input.code:q} --samples {input.samples:q} --mapping {input.mapping:q} "
        "--protein-dir {params.proteins:q} --outdir {output.root:q} "
        "--batch-size {params.batch_size} > {log:q} 2>&1"


rule annotate_kofam:
    input:
        protein=kegg_batch_input,
        reference=f"{KEGG_REFERENCE}/reference.json",
        verified=f"{KEGG}/reference_qc.json",
        code=f"{SCRIPTS}/run_kofam.py",
        cleanup=f"{SCRIPTS}/cleanup_work.py",
        verifier=f"{SCRIPTS}/verify_kegg_reference.py",
        preparation=f"{SCRIPTS}/prepare_kegg_reference.py",
        common=f"{SCRIPTS}/common.py"
    output:
        detail=f"{KEGG_ANNOTATIONS}/{{batch}}/detail.tsv",
        hits=f"{KEGG_ANNOTATIONS}/{{batch}}/gene_kos.tsv",
        genes=f"{KEGG_ANNOTATIONS}/{{batch}}/genes.tsv",
        provenance=f"{KEGG_ANNOTATIONS}/{{batch}}/provenance.json"
    params:
        outdir=lambda wc: f"{KEGG_ANNOTATIONS}/{wc.batch}",
        workdir=lambda wc: f"{WORK}/kegg/{wc.batch}",
        command="exec_annotation",
        keep=["--keep-intermediates"] if config["storage"]["keep_intermediates"] else []
    threads: 4
    resources: mem_mb=8000
    log: f"{LOG}/kegg/annotation/{{batch}}.log"
    benchmark: f"{LOG}/kegg/benchmarks/{{batch}}.tsv"
    conda: "../envs/kofam.yaml"
    shell:
        "{PYTHON:q} {input.code:q} --protein {input.protein:q} --species orthogroup_representatives "
        "--reference {input.reference:q} --output-dir {params.outdir:q} --work-dir {params.workdir:q} "
        "--command {params.command:q} --threads {threads} {params.keep:q} > {log:q} 2>&1"


rule assign_og_kos:
    input:
        representatives=lambda wc: checkpoints.select_kegg_representatives.get().output.root,
        genes=lambda wc: kegg_batch_outputs(wc, "genes.tsv"),
        hits=lambda wc: kegg_batch_outputs(wc, "gene_kos.tsv"),
        details=lambda wc: kegg_batch_outputs(wc, "detail.tsv"),
        provenance=lambda wc: kegg_batch_outputs(wc, "provenance.json"),
        verified=f"{KEGG}/reference_qc.json",
        code=f"{SCRIPTS}/og_kegg.py",
        helpers=[f"{SCRIPTS}/{name}.py" for name in [
            "select_kegg_representatives", "aggregate_ko_tpm", "common"]]
    output:
        hits=f"{KEGG}/og_kos.tsv",
        groups=f"{KEGG}/orthogroups.tsv",
        provenance=f"{KEGG}/annotation_provenance.json"
    params: annotation=KEGG_ANNOTATIONS, outdir=KEGG
    log: f"{LOG}/kegg/assign.log"
    conda: "../envs/analysis.yaml"
    resources: mem_mb=8000
    shell:
        "{PYTHON:q} {input.code:q} assign --representatives {input.representatives:q} "
        "--annotation-dir {params.annotation:q} --outdir {params.outdir:q} > {log:q} 2>&1"


rule merge_kegg:
    input:
        samples=f"{META}/samples.tsv",
        groups=f"{KEGG}/orthogroups.tsv",
        hits=f"{KEGG}/og_kos.tsv",
        provenance=f"{KEGG}/annotation_provenance.json",
        tables=lambda wc: run_outputs(wc, "tsv"),
        reports=lambda wc: run_outputs(wc, "qc.json"),
        code=f"{SCRIPTS}/og_kegg.py",
        helpers=[f"{SCRIPTS}/{name}.py" for name in [
            "select_kegg_representatives", "aggregate_ko_tpm", "common"]]
    output:
        summed=f"{KEGG}/ko_tpm_sum.tsv",
        wide=f"{KEGG}/ko_tpm_sum_wide.tsv",
        support=f"{KEGG}/ko_support.tsv",
        qc=f"{KEGG}/mapping_qc.tsv",
        provenance=f"{KEGG}/provenance.json"
    params:
        annotation=KEGG, runs=f"{TPM}/runs", outdir=KEGG,
        ambiguity=config["kegg"]["ambiguity"]
    log: f"{LOG}/kegg/merge.log"
    benchmark: f"{LOG}/kegg/benchmarks/merge.tsv"
    conda: "../envs/analysis.yaml"
    resources: mem_mb=8000
    shell:
        "{PYTHON:q} {input.code:q} aggregate --samples {input.samples:q} "
        "--annotation-dir {params.annotation:q} --og-run-dir {params.runs:q} "
        "--outdir {params.outdir:q} --ambiguity {params.ambiguity:q} > {log:q} 2>&1"
