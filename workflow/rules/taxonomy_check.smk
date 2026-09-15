rule check_taxonomy:
    wildcard_constraints: phylo_branch=MOLECULAR_BRANCH_PATTERN
    input:
        tree=f"{PHYLO_RUN}/species_tree.nwk",
        tree_qc=f"{PHYLO_RUN}/species_tree.json",
        samples=lambda wc: f"{PHENOTYPED}/selection/samples.tsv" if wc.phylo_branch == PHYLO_BRANCHES["phenotyped"] else f"{META}/samples.tsv",
        taxonomy=TAXONOMY_DB,
        code=f"{SCRIPTS}/taxonomy_check.py",
        monophy=f"{SCRIPTS}/taxonomy_check.R",
        common=f"{SCRIPTS}/common.py"
    output: report=directory(f"{PHYLO_RUN}/taxonomy_check")
    params: settings=json.dumps(TAXONOMY_CHECK, sort_keys=True)
    conda: "../envs/monophy.yaml"
    resources: mem_mb=8000
    log: f"{LOG}/{{phylo_branch}}/taxonomy_check.log"
    shell:
        "{PYTHON:q} {input.code:q} --tree {input.tree:q} --tree-qc {input.tree_qc:q} "
        "--samples {input.samples:q} --taxonomy {input.taxonomy:q} "
        "--outdir {output.report:q} --settings {params.settings:q} > {log:q} 2>&1"
