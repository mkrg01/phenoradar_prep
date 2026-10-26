if config.get("build_manifest"):
    rule import_build_expression:
        input:
            completion=config['build_manifest'],
            samples=f"{META}/samples.tsv",
            expression=lambda wc: BUILD_INPUTS.product(wc.run, 'run')['expression']['path'],
            qc=lambda wc: BUILD_INPUTS.product(wc.run, 'run')['expression_qc']['path'],
            code=ancient(f"{SCRIPTS}/build_products.py"),
            helpers=[f"{SCRIPTS}/portable_build.py", f"{SCRIPTS}/dataset_assets.py", f"{SCRIPTS}/product_storage.py"]
        output:
            tpm=f"{TPM}/runs/{{run}}.tsv.gz",
            qc=f"{TPM}/runs/{{run}}.qc.json"
        conda: "../envs/analysis.yaml"
        params: index=BUILD_INDEX, implementation=BUILD_IMPORT_SIGNATURE
        resources: mem_mb=4000
        shell:
            "{PYTHON:q} {input.code:q} expression --completion {input.completion:q} "
            "--samples {input.samples:q} --run {wildcards.run:q} --output {output.tpm:q} "
            "--qc {output.qc:q} --index {params.index:q}"
else:
    from sample_table import expression_batches

    def expression_plan(wc):
        return expression_batches(checkpoints.select_metadata.get().output.samples)

    def expression_batch(wc):
        rows = expression_plan(wc)[0].get(wc.batch)
        if rows is None:
            raise WorkflowError(f'unknown expression batch: {wc.batch}')
        return rows

    rule aggregate_tpm:
        input:
            samples=f"{META}/samples.tsv",
            abundance=lambda wc: [row['abundance'] for row in expression_batch(wc)],
            mapping=f"{MAPPING}/snapshot.json",
            tables=f"{MAPPING}/species",
            code=f"{SCRIPTS}/aggregate_tpm.py",
            common=[f"{SCRIPTS}/{name}.py" for name in
                    ("common", "mapping_tables", "dataset_assets", "product_storage", "sample_table")]
        output: directory(f"{TPM}/batches/{{batch}}")
        wildcard_constraints: batch=r"batch_[0-9]+"
        params:
            cache=config.get("expression_cache") or "",
            runs=lambda wc: [row['run'] for row in expression_batch(wc)]
        log: f"{LOG}/{ORTHOGROUP_EXPRESSION}/{{batch}}.log"
        conda: "../envs/analysis.yaml"
        resources: mem_mb=2000
        shell:
            "{PYTHON:q} {input.code:q} --samples {input.samples:q} --runs {params.runs:q} "
            "--mapping {input.mapping:q} --outdir {output:q} "
            "--cache-dir={params.cache:q} --gzip > {log:q} 2>&1"

    localrules: publish_tpm_run

    rule publish_tpm_run:
        input:
            batch=lambda wc: f'{TPM}/batches/{expression_plan(wc)[1][wc.run]}',
            code=f'{SCRIPTS}/aggregate_tpm.py',
            helpers=[f'{SCRIPTS}/{name}.py' for name in
                     ('common', 'mapping_tables', 'dataset_assets', 'product_storage', 'sample_table')]
        output:
            tpm=f'{TPM}/runs/{{run}}.tsv.gz',
            qc=f'{TPM}/runs/{{run}}.qc.json'
        conda: '../envs/analysis.yaml'
        resources: mem_mb=256
        shell:
            '{PYTHON:q} {input.code:q} --batch-dir {input.batch:q} --run {wildcards.run:q} '
            '--output {output.tpm:q} --qc {output.qc:q}'


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
    resources: mem_mb=32000
    shell:
        "{PYTHON:q} {input.code:q} --samples {input.samples:q} --run-dir {params.runs:q} "
        "--outdir {params.out:q} > {log:q} 2>&1"
