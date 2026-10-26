from alignment_batches import batch_targets


def alignment_targets(wc):
    plan = checkpoints.plan_alignment_batches.get().output.plan
    return ALIGNMENT_FINAL + batch_targets(plan, ALIGNMENTS, f'{LOG}/{ORTHOGROUP_ALIGNMENTS}')


def alignment_batch_file(wc):
    plan = checkpoints.plan_alignment_batches.get().output.plan
    return str(Path(plan) / f'{wc.batch}.json')


checkpoint collect_orthogroup_proteins:
    input:
        samples=f"{META}/samples.tsv",
        mapping=f"{MAPPING}/snapshot.json",
        tables=f"{MAPPING}/species",
        proteins=lambda wc: sorted({f'{PROTEINS}/{r["odb_species"]}_protein.fa.gz' for r in sample_rows(wc)}),
        code=f"{SCRIPTS}/align_orthogroups.py",
        helpers=[f"{SCRIPTS}/common.py", f"{SCRIPTS}/busco_phylogeny.py"]
    output: fasta=directory(ALIGNMENT_INPUTS)
    params: proteins=PROTEINS
    resources: mem_mb=8000
    conda: "../envs/alignment.yaml"
    log: f"{LOG}/{ORTHOGROUP_ALIGNMENTS}/collect.log"
    shell:
        "{PYTHON:q} {input.code:q} collect --samples {input.samples:q} --mapping {input.mapping:q} "
        "--protein-dir {params.proteins:q} --outdir {output.fasta:q} > {log:q} 2>&1"


checkpoint plan_alignment_batches:
    input:
        proteins=lambda wc: checkpoints.collect_orthogroup_proteins.get().output.fasta,
        code=f'{SCRIPTS}/alignment_batches.py'
    output: plan=directory(f'{WORK}/{ORTHOGROUP_ALIGNMENTS}/batches')
    resources: mem_mb=4000
    conda: '../envs/alignment.yaml'
    shell:
        '{PYTHON:q} {input.code:q} plan --inputs {input.proteins:q} --outdir {output.plan:q}'


rule align_orthogroup_batch:
    wildcard_constraints: batch='batch_[0-9]+'
    input:
        batch=alignment_batch_file,
        code=f'{SCRIPTS}/alignment_batches.py',
        aligner=f'{SCRIPTS}/align_orthogroups.py'
    output: receipt=f'{ALIGNMENTS}/batches/{{batch}}.json'
    params: inputs=ALIGNMENT_INPUTS, outdir=ALIGNMENTS, reports=f'{LOG}/{ORTHOGROUP_ALIGNMENTS}'
    threads: 12
    resources: mem_mb=32000
    conda: '../envs/alignment.yaml'
    log: f'{LOG}/{ORTHOGROUP_ALIGNMENTS}/batches/{{batch}}.log'
    benchmark: f'{LOG}/{ORTHOGROUP_ALIGNMENTS}/batches/benchmarks/{{batch}}.tsv'
    shell:
        '{PYTHON:q} {input.code:q} batch --batch {input.batch:q} --inputs {params.inputs:q} '
        '--outdir {params.outdir:q} --reports {params.reports:q} --receipt {output.receipt:q} '
        '--cpus {threads} > {log:q} 2>&1'


rule finish_alignments:
    input:
        proteins=lambda wc: checkpoints.collect_orthogroup_proteins.get().output.fasta,
        batches=lambda wc: batch_targets(checkpoints.plan_alignment_batches.get().output.plan,
                                        ALIGNMENTS, f'{LOG}/{ORTHOGROUP_ALIGNMENTS}'),
        code=f"{SCRIPTS}/align_orthogroups.py",
        helpers=[f"{SCRIPTS}/common.py", f"{SCRIPTS}/busco_phylogeny.py"]
    output: provenance=f"{ALIGNMENTS}/provenance.json"
    params: outdir=ALIGNMENTS, reports=f"{LOG}/{ORTHOGROUP_ALIGNMENTS}"
    resources: mem_mb=8000
    conda: "../envs/alignment.yaml"
    log: f"{LOG}/{ORTHOGROUP_ALIGNMENTS}/finish.log"
    shell:
        "{PYTHON:q} {input.code:q} finish --inputs {input.proteins:q} --outdir {params.outdir:q} "
        "--reports {params.reports:q} > {log:q} 2>&1"
