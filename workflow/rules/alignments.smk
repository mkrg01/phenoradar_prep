def alignment_files(wc, suffix, root):
    inputs = checkpoints.collect_orthogroup_proteins.get().output.fasta
    with open(Path(inputs) / "provenance.json") as handle:
        groups = json.load(handle)["orthogroups"]
    return [f"{root}/{og}.{suffix}" for og in groups]


def alignment_targets(wc):
    # Request each OG directly: Snakemake otherwise considers a deleted MSA an
    # unnecessary intermediate when the final provenance still exists.
    return (ALIGNMENT_FINAL + alignment_files(wc, "faa", ALIGNMENTS)
            + alignment_files(wc, "json", f"{LOG}/alignments"))


checkpoint collect_orthogroup_proteins:
    input:
        samples=f"{META}/samples.tsv",
        database=f"{MERGED}/mappings.sqlite",
        proteins=lambda wc: sorted({f'{PROTEINS}/{r["odb_species"]}_protein.fa' for r in sample_rows(wc)}),
        code=f"{SCRIPTS}/align_orthogroups.py",
        helpers=[f"{SCRIPTS}/common.py", f"{SCRIPTS}/busco_phylogeny.py"]
    output: fasta=directory(ALIGNMENT_INPUTS)
    params: proteins=PROTEINS
    resources: mem_mb=config["alignment"]["mem_gb"] * 1000
    conda: "../envs/alignment.yaml"
    log: f"{LOG}/alignments/collect.log"
    shell:
        "{PYTHON:q} {input.code:q} collect --samples {input.samples:q} --database {input.database:q} "
        "--protein-dir {params.proteins:q} --outdir {output.fasta:q} > {log:q} 2>&1"


rule align_orthogroup:
    wildcard_constraints: og="[A-Za-z0-9][A-Za-z0-9_.-]*"
    input:
        proteins=lambda wc: checkpoints.collect_orthogroup_proteins.get().output.fasta,
        code=f"{SCRIPTS}/align_orthogroups.py",
        helpers=[f"{SCRIPTS}/common.py", f"{SCRIPTS}/busco_phylogeny.py"]
    output:
        alignment=f"{ALIGNMENTS}/{{og}}.faa",
        provenance=f"{LOG}/alignments/{{og}}.json"
    params: fasta=lambda wc: f"{ALIGNMENT_INPUTS}/{wc.og}.faa"
    threads: config["alignment"]["threads"]
    resources: mem_mb=config["alignment"]["mem_gb"] * 1000
    conda: "../envs/alignment.yaml"
    log: f"{LOG}/alignments/{{og}}.log"
    benchmark: f"{LOG}/benchmarks/alignment_{{og}}.tsv"
    shell:
        "{PYTHON:q} {input.code:q} align --fasta {params.fasta:q} --output {output.alignment:q} "
        "--provenance {output.provenance:q} --threads {threads} > {log:q} 2>&1"


rule finish_alignments:
    input:
        proteins=lambda wc: checkpoints.collect_orthogroup_proteins.get().output.fasta,
        alignments=lambda wc: alignment_files(wc, "faa", ALIGNMENTS),
        reports=lambda wc: alignment_files(wc, "json", f"{LOG}/alignments"),
        code=f"{SCRIPTS}/align_orthogroups.py",
        helpers=[f"{SCRIPTS}/common.py", f"{SCRIPTS}/busco_phylogeny.py"]
    output:
        provenance=f"{ALIGNMENTS}/provenance.json"
    params: outdir=ALIGNMENTS, reports=f"{LOG}/alignments"
    resources: mem_mb=config["alignment"]["mem_gb"] * 1000
    conda: "../envs/alignment.yaml"
    log: f"{LOG}/alignments/finish.log"
    shell:
        "{PYTHON:q} {input.code:q} finish --inputs {input.proteins:q} --outdir {params.outdir:q} "
        "--reports {params.reports:q} > {log:q} 2>&1"
