rule kegg_verify_reference:
    input:
        reference=f"{KEGG_REFERENCE}/reference.json",
        inventory=f"{KEGG_REFERENCE}/files.json",
        modules=f"{KEGG_REFERENCE}/ko_modules.tsv",
        pathways=f"{KEGG_REFERENCE}/ko_pathways.tsv",
        code=f"{SCRIPTS}/publish_kegg_reference.py",
        verifier=f"{SCRIPTS}/verify_kegg_reference.py",
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


rule annotate_kofam:
    input:
        protein=f"{PROTEINS}/{{species}}_protein.fa",
        reference=f"{KEGG_REFERENCE}/reference.json",
        verified=f"{KEGG}/reference_qc.json",
        code=f"{SCRIPTS}/run_kofam.py",
        verifier=f"{SCRIPTS}/verify_kegg_reference.py",
        preparation=f"{SCRIPTS}/prepare_kegg_reference.py",
        translation=f"{SCRIPTS}/translate_cds.py",
        common=f"{SCRIPTS}/common.py"
    output:
        detail=f"{KEGG_SPECIES}/{{species}}/detail.tsv",
        hits=f"{KEGG_SPECIES}/{{species}}/gene_kos.tsv",
        genes=f"{KEGG_SPECIES}/{{species}}/genes.tsv",
        provenance=f"{KEGG_SPECIES}/{{species}}/provenance.json"
    params:
        species=lambda wc: species_row(wc)["species"],
        outdir=lambda wc: f"{KEGG_SPECIES}/{wc.species}",
        workdir=lambda wc: f"{WORK}/kegg/{wc.species}",
        command=config["kegg"]["command"]
    threads: config["kegg"]["threads"]
    resources: mem_mb=config["kegg"]["mem_gb"] * 1000
    log: f"{LOG}/kegg/species/{{species}}.log"
    benchmark: f"{KEGG_SPECIES}/{{species}}/benchmark.tsv"
    conda: "../envs/kofam.yaml"
    shell:
        "{PYTHON:q} {input.code:q} --protein {input.protein:q} --species {params.species:q} "
        "--reference {input.reference:q} --output-dir {params.outdir:q} --work-dir {params.workdir:q} "
        "--command {params.command:q} --threads {threads} > {log:q} 2>&1"


rule aggregate_ko_tpm:
    input:
        samples=f"{META}/samples.tsv",
        abundance=lambda wc: run_row(wc)["abundance"],
        genes=lambda wc: f'{KEGG_SPECIES}/{run_row(wc)["odb_species"]}/genes.tsv',
        hits=lambda wc: f'{KEGG_SPECIES}/{run_row(wc)["odb_species"]}/gene_kos.tsv',
        detail=lambda wc: f'{KEGG_SPECIES}/{run_row(wc)["odb_species"]}/detail.tsv',
        provenance=lambda wc: f'{KEGG_SPECIES}/{run_row(wc)["odb_species"]}/provenance.json',
        code=f"{SCRIPTS}/aggregate_ko_tpm.py",
        common=f"{SCRIPTS}/common.py"
    output:
        tpm=f"{KEGG_RUNS}/{{run}}.tsv",
        qc=f"{KEGG_RUNS}/{{run}}.qc.json"
    params:
        annotation=lambda wc: f'{KEGG_SPECIES}/{run_row(wc)["odb_species"]}',
        ambiguity=config["kegg"]["ambiguity"]
    log: f"{LOG}/kegg/tpm/{{run}}.log"
    conda: "../envs/analysis.yaml"
    resources: mem_mb=2000
    shell:
        "{PYTHON:q} {input.code:q} --samples {input.samples:q} --run {wildcards.run:q} "
        "--annotation-dir {params.annotation:q} --output {output.tpm:q} --qc {output.qc:q} "
        "--ambiguity {params.ambiguity:q} > {log:q} 2>&1"


rule merge_kegg:
    input:
        samples=f"{META}/samples.tsv",
        genes=lambda wc: kegg_species_outputs(wc, "genes.tsv"),
        hits=lambda wc: kegg_species_outputs(wc, "gene_kos.tsv"),
        details=lambda wc: kegg_species_outputs(wc, "detail.tsv"),
        provenance=lambda wc: kegg_species_outputs(wc, "provenance.json"),
        tables=lambda wc: kegg_run_outputs(wc, "tsv"),
        reports=lambda wc: kegg_run_outputs(wc, "qc.json"),
        code=f"{SCRIPTS}/merge_kegg.py",
        aggregation=f"{SCRIPTS}/aggregate_ko_tpm.py",
        common=f"{SCRIPTS}/common.py"
    output:
        hits=f"{KEGG}/gene_kos.tsv",
        genes=f"{KEGG}/genes.tsv",
        summed=f"{KEGG}/ko_tpm_sum.tsv",
        wide=f"{KEGG}/ko_tpm_sum_wide.tsv",
        support=f"{KEGG}/ko_support.tsv",
        qc=f"{KEGG}/mapping_qc.tsv"
    params: annotation=KEGG_SPECIES, runs=KEGG_RUNS, outdir=KEGG
    log: f"{LOG}/kegg/merge.log"
    conda: "../envs/analysis.yaml"
    resources: mem_mb=8000
    shell:
        "{PYTHON:q} {input.code:q} --samples {input.samples:q} --annotation-dir {params.annotation:q} "
        "--run-dir {params.runs:q} --outdir {params.outdir:q} > {log:q} 2>&1"
