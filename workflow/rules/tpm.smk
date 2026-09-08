rule aggregate_tpm:
    input:
        samples=f"{META}/samples.tsv",
        abundance=lambda wc: run_row(wc)["abundance"],
        database=f"{MERGED}/mappings.sqlite",
        code=f"{SCRIPTS}/aggregate_tpm.py",
        common=f"{SCRIPTS}/common.py"
    output:
        tpm=f"{TPM}/runs/{{run}}.tsv",
        qc=f"{TPM}/runs/{{run}}.qc.json"
    params: multimap=config["tpm"]["multimap"]
    log: f"{LOG}/tpm/{{run}}.log"
    conda: "../envs/analysis.yaml"
    resources: mem_mb=2000
    shell:
        "{PYTHON:q} {input.code:q} --samples {input.samples:q} --run {wildcards.run:q} "
        "--database {input.database:q} --output {output.tpm:q} --qc {output.qc:q} "
        "--multimap {params.multimap:q} > {log:q} 2>&1"


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
    log: f"{LOG}/merge_tpm.log"
    conda: "../envs/analysis.yaml"
    resources: mem_mb=8000
    shell:
        "{PYTHON:q} {input.code:q} --samples {input.samples:q} --run-dir {params.runs:q} "
        "--outdir {params.out:q} > {log:q} 2>&1"
