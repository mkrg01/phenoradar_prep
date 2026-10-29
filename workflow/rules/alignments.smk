from alignment_batches import batch_targets

SELECTED_ALIGNMENT_INPUTS = f'{WORK}/{ORTHOGROUP_ALIGNMENTS}/selected'
ALIGNMENT_OCCUPANCY = f'{WORK}/{ORTHOGROUP_ALIGNMENTS}/occupancy'
ALIGNMENT_CACHE = f'{WORK}/{ORTHOGROUP_ALIGNMENTS}/cache'


def alignment_targets(wc):
    plan = checkpoints.plan_alignment_batches.get().output.plan
    return (ALIGNMENT_FINAL + [f'{ALIGNMENTS}/occupancy.tsv']
            + batch_targets(plan, ALIGNMENTS, f'{LOG}/{ORTHOGROUP_ALIGNMENTS}'))


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


rule measure_alignment_occupancy:
    input:
        proteins=lambda wc: checkpoints.collect_orthogroup_proteins.get().output.fasta,
        samples=f'{META}/samples.tsv', mapping=f'{MAPPING}/snapshot.json',
        tables=f'{MAPPING}/species', code=f'{SCRIPTS}/alignment_selection.py'
    output:
        counts=f'{ALIGNMENT_OCCUPANCY}/occupancy.tsv',
        provenance=f'{ALIGNMENT_OCCUPANCY}/provenance.json'
    params: outdir=ALIGNMENT_OCCUPANCY
    resources: mem_mb=4000
    conda: '../envs/alignment.yaml'
    log: f'{LOG}/{ORTHOGROUP_ALIGNMENTS}/occupancy.log'
    shell:
        '{PYTHON:q} {input.code:q} measure --inputs {input.proteins:q} --samples {input.samples:q} '
        '--mapping {input.mapping:q} --outdir {params.outdir:q} > {log:q} 2>&1'


checkpoint plan_alignment_batches:
    input:
        proteins=lambda wc: checkpoints.collect_orthogroup_proteins.get().output.fasta,
        counts=rules.measure_alignment_occupancy.output.counts,
        occupancy=rules.measure_alignment_occupancy.output.provenance,
        code=f'{SCRIPTS}/alignment_batches.py', selection_code=f'{SCRIPTS}/alignment_selection.py'
    output:
        fasta=directory(SELECTED_ALIGNMENT_INPUTS),
        plan=directory(f'{WORK}/{ORTHOGROUP_ALIGNMENTS}/batches')
    params: occupancy=ALIGNMENT_OCCUPANCY, fraction=config['alignment']['min_sample_fraction']
    resources: mem_mb=4000
    conda: '../envs/alignment.yaml'
    log: f'{LOG}/{ORTHOGROUP_ALIGNMENTS}/selection.log'
    shell:
        '({PYTHON:q} {input.selection_code:q} select --inputs {input.proteins:q} '
        '--occupancy {params.occupancy:q} --outdir {output.fasta:q} '
        '--min-sample-fraction {params.fraction} && '
        '{PYTHON:q} {input.code:q} plan --inputs {output.fasta:q} --outdir {output.plan:q}) > {log:q} 2>&1'


rule align_orthogroup_batch:
    wildcard_constraints: batch='batch_[0-9]+'
    input:
        batch=alignment_batch_file,
        selection_manifest=f'{SELECTED_ALIGNMENT_INPUTS}/provenance.json',
        code=f'{SCRIPTS}/alignment_batches.py',
        aligner=f'{SCRIPTS}/align_orthogroups.py'
    output: receipt=f'{ALIGNMENTS}/batches/{{batch}}.json'
    params: inputs=SELECTED_ALIGNMENT_INPUTS, outdir=ALIGNMENTS, reports=f'{LOG}/{ORTHOGROUP_ALIGNMENTS}', cache=ALIGNMENT_CACHE
    threads: 12
    resources: mem_mb=32000
    conda: '../envs/alignment.yaml'
    log: f'{LOG}/{ORTHOGROUP_ALIGNMENTS}/batches/{{batch}}.log'
    benchmark: f'{LOG}/{ORTHOGROUP_ALIGNMENTS}/batches/benchmarks/{{batch}}.tsv'
    shell:
        '{PYTHON:q} {input.code:q} batch --batch {input.batch:q} --inputs {params.inputs:q} '
        '--outdir {params.outdir:q} --reports {params.reports:q} --receipt {output.receipt:q} '
        '--cpus {threads} --cache {params.cache:q} > {log:q} 2>&1'


rule finish_alignments:
    input:
        proteins=lambda wc: checkpoints.plan_alignment_batches.get().output.fasta,
        batches=lambda wc: batch_targets(checkpoints.plan_alignment_batches.get().output.plan,
                                        ALIGNMENTS, f'{LOG}/{ORTHOGROUP_ALIGNMENTS}'),
        code=f"{SCRIPTS}/alignment_selection.py",
        helpers=[f"{SCRIPTS}/common.py", f"{SCRIPTS}/busco_phylogeny.py", f"{SCRIPTS}/alignment_batches.py"]
    output:
        provenance=f"{ALIGNMENTS}/provenance.json",
        occupancy=f"{ALIGNMENTS}/occupancy.tsv"
    params: outdir=ALIGNMENTS, reports=f"{LOG}/{ORTHOGROUP_ALIGNMENTS}", cache=ALIGNMENT_CACHE
    resources: mem_mb=8000
    conda: "../envs/alignment.yaml"
    log: f"{LOG}/{ORTHOGROUP_ALIGNMENTS}/finish.log"
    shell:
        "{PYTHON:q} {input.code:q} finish --inputs {input.proteins:q} --outdir {params.outdir:q} "
        "--reports {params.reports:q} --cache {params.cache:q} > {log:q} 2>&1"
