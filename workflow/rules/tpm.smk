if config.get("build_manifest") and COMPLETED_BUILD['schema_version'] >= 5:
    if config['tpm'] != COMPLETED_BUILD['tpm']:
        raise WorkflowError('downstream TPM settings differ from the database')
    BUILD_RUNS = {p['row']['run']: p for p in COMPLETED_BUILD['products'].values()}

    rule import_build_expression:
        input:
            completion=config['build_manifest'],
            samples=f"{META}/samples.tsv",
            expression=lambda wc: BUILD_RUNS[wc.run]['expression']['path'],
            qc=lambda wc: BUILD_RUNS[wc.run]['expression_qc']['path'],
            code=f"{SCRIPTS}/build_products.py",
            helpers=[f"{SCRIPTS}/portable_build.py", f"{SCRIPTS}/dataset_assets.py"]
        output:
            tpm=f"{TPM}/runs/{{run}}.tsv",
            qc=f"{TPM}/runs/{{run}}.qc.json"
        params: multimap=config['tpm']['multimap']
        conda: "../envs/analysis.yaml"
        resources: mem_mb=2000
        shell:
            "{PYTHON:q} {input.code:q} expression --completion {input.completion:q} "
            "--samples {input.samples:q} --run {wildcards.run:q} --output {output.tpm:q} "
            "--qc {output.qc:q} --multimap {params.multimap:q}"
else:
    rule aggregate_tpm:
        input:
            samples=f"{META}/samples.tsv",
            abundance=lambda wc: run_row(wc)["abundance"],
            mapping=f"{MAPPING}/snapshot.json",
            tables=f"{MAPPING}/species",
            code=f"{SCRIPTS}/aggregate_tpm.py",
            common=[f"{SCRIPTS}/common.py", f"{SCRIPTS}/mapping_tables.py", f"{SCRIPTS}/dataset_assets.py"]
        output:
            tpm=f"{TPM}/runs/{{run}}.tsv",
            qc=f"{TPM}/runs/{{run}}.qc.json"
        params:
            multimap=config["tpm"]["multimap"],
            cache=config.get("expression_cache") or ""
        log: f"{LOG}/{ORTHOGROUP_EXPRESSION}/{{run}}.log"
        conda: "../envs/analysis.yaml"
        resources: mem_mb=2000
        shell:
            "{PYTHON:q} {input.code:q} --samples {input.samples:q} --run {wildcards.run:q} "
            "--mapping {input.mapping:q} --output {output.tpm:q} --qc {output.qc:q} "
            "--multimap {params.multimap:q} --cache-dir={params.cache:q} > {log:q} 2>&1"


rule merge_tpm:
    input:
        samples=f"{META}/samples.tsv",
        tables=lambda wc: run_outputs(wc, "tsv"),
        reports=lambda wc: run_outputs(wc, "qc.json"),
        code=f"{SCRIPTS}/merge_tpm.py",
        common=f"{SCRIPTS}/common.py"
    output:
        normalized=f"{TPM}/tpm.tsv",
        normalized_wide=f"{TPM}/tpm_wide.tsv",
        summed=f"{TPM}/tpm_sum.tsv",
        summed_wide=f"{TPM}/tpm_sum_wide.tsv",
        qc=f"{TPM}/mapping_qc.tsv"
    params: runs=f"{TPM}/runs", out=TPM
    log: f"{LOG}/{ORTHOGROUP_EXPRESSION}/merge.log"
    conda: "../envs/analysis.yaml"
    resources: mem_mb=8000
    shell:
        "{PYTHON:q} {input.code:q} --samples {input.samples:q} --run-dir {params.runs:q} "
        "--outdir {params.out:q} > {log:q} 2>&1"
