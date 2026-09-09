def phylogeny_tables(wc):
    if not PHY["busco_full_dir"]:
        raise WorkflowError("set phylogeny.busco_full_dir to per-species BUSCO full tables; see docs/phylogeny.md")
    if not PHY["outgroup"]:
        raise WorkflowError("set phylogeny.outgroup to an exact selected species label")
    return sorted({str(Path(PHY["busco_full_dir"]) / (r["species"] + PHY["busco_full_suffix"])) for r in sample_rows(wc)})


def phylogeny_plan_rows(wc, table):
    output = checkpoints.plan_phylogeny.get().output
    with open(getattr(output, table)) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def phylogeny_species_row(wc):
    rows = [r for r in phylogeny_plan_rows(wc, "species") if r["species"] == wc.species]
    if len(rows) != 1:
        raise WorkflowError(f"species absent from phylogeny plan: {wc.species}")
    return rows[0]


def phylogeny_marker_files(wc, suffix):
    return [f'{PHYLO}/gene_trees/{r["marker"]}.{suffix}' for r in phylogeny_plan_rows(wc, "markers")]


def dating_calibrations(wc):
    if PHY["dating"]["calibration_source"] == "timetree":
        return f"{PHYLO}/timetree/calibrations.tsv"
    if not PHY["dating"]["calibrations"]:
        raise WorkflowError("absolute dating requires phylogeny.dating.calibrations; no arbitrary root age is used")
    return PHY["dating"]["calibrations"]


checkpoint plan_phylogeny:
    input:
        samples=f"{META}/samples.tsv",
        tables=phylogeny_tables,
        code=f"{SCRIPTS}/busco_phylogeny.py",
        common=f"{SCRIPTS}/common.py"
    output:
        species=f"{PHYLO_PLAN}/species.tsv",
        markers=f"{PHYLO_PLAN}/markers.tsv",
        stats=f"{PHYLO_PLAN}/marker_stats.tsv",
        provenance=f"{PHYLO_PLAN}/provenance.json"
    params:
        outdir=PHYLO_PLAN,
        settings=json.dumps({k: PHY[k] for k in ["outgroup", "busco_full_dir", "busco_full_suffix", "lineage",
            "sequence_dir", "sequence_suffix", "min_taxa", "max_markers"]}, sort_keys=True)
    log: f"{LOG}/phylogeny/plan.log"
    conda: "../envs/phylogeny.yaml"
    resources: mem_mb=PHY["preparation_mem_gb"] * 1000
    shell:
        "{PYTHON:q} {input.code:q} plan --samples {input.samples:q} "
        "--outdir {params.outdir:q} --settings {params.settings:q} > {log:q} 2>&1"


rule extract_busco_proteins:
    input:
        markers=f"{PHYLO_PLAN}/markers.tsv",
        table=lambda wc: phylogeny_species_row(wc)["busco_table"],
        sequences=lambda wc: phylogeny_species_row(wc)["sequences"],
        code=f"{SCRIPTS}/busco_phylogeny.py",
        common=f"{SCRIPTS}/common.py"
    output:
        proteins=f"{PHYLO}/species/{{species}}.faa",
        qc=f"{PHYLO}/species/{{species}}.json"
    params:
        settings=json.dumps({k: PHY[k] for k in ["lineage", "sequence_mode", "translation_table",
            "min_protein_length", "max_unknown_fraction"]}, sort_keys=True)
    conda: "../envs/phylogeny.yaml"
    resources: mem_mb=PHY["preparation_mem_gb"] * 1000
    log: f"{LOG}/phylogeny/extract/{{species}}.log"
    shell:
        "{PYTHON:q} {input.code:q} extract --species {wildcards.species:q} --table {input.table:q} "
        "--sequences {input.sequences:q} --markers {input.markers:q} --output {output.proteins:q} "
        "--qc {output.qc:q} --settings {params.settings:q} > {log:q} 2>&1"


rule collect_busco_markers:
    input:
        manifest=f"{PHYLO_PLAN}/species.tsv",
        markers=f"{PHYLO_PLAN}/markers.tsv",
        proteins=lambda wc: [f'{PHYLO}/species/{r["species"]}.faa' for r in phylogeny_plan_rows(wc, "species")],
        reports=lambda wc: [f'{PHYLO}/species/{r["species"]}.json' for r in phylogeny_plan_rows(wc, "species")],
        code=f"{SCRIPTS}/busco_phylogeny.py",
        common=f"{SCRIPTS}/common.py"
    output: fasta=directory(f"{PHYLO}/markers")
    params: species_dir=f"{PHYLO}/species"
    conda: "../envs/phylogeny.yaml"
    resources: mem_mb=PHY["preparation_mem_gb"] * 1000
    log: f"{LOG}/phylogeny/collect.log"
    shell:
        "{PYTHON:q} {input.code:q} collect --manifest {input.manifest:q} --markers {input.markers:q} "
        "--species-dir {params.species_dir:q} --outdir {output.fasta:q} > {log:q} 2>&1"


rule align_busco_marker:
    wildcard_constraints: marker="[A-Za-z0-9][A-Za-z0-9_.-]*"
    input:
        markers=rules.collect_busco_markers.output.fasta,
        code=f"{SCRIPTS}/infer_phylogeny.py",
        helpers=[f"{SCRIPTS}/busco_phylogeny.py", f"{SCRIPTS}/common.py"]
    output:
        alignment=f"{PHYLO}/alignments/raw/{{marker}}.faa",
        qc=f"{PHYLO}/alignments/raw/{{marker}}.json"
    params:
        fasta=lambda wc: f"{PHYLO}/markers/{wc.marker}.faa",
        command=PHY["famsa_command"],
        settings=json.dumps({"min_taxa": PHY["min_taxa"]}, sort_keys=True)
    threads: PHY["align_threads"]
    resources: mem_mb=PHY["alignment_mem_gb"] * 1000
    conda: "../envs/phylogeny.yaml"
    log: f"{LOG}/phylogeny/align/{{marker}}.log"
    benchmark: f"{PHYLO}/benchmarks/align.{{marker}}.tsv"
    shell:
        "{PYTHON:q} {input.code:q} align --fasta {params.fasta:q} --output {output.alignment:q} "
        "--qc {output.qc:q} --command {params.command:q} --threads {threads} "
        "--settings {params.settings:q} > {log:q} 2>&1"


rule trim_busco_marker:
    wildcard_constraints: marker="[A-Za-z0-9][A-Za-z0-9_.-]*"
    input:
        alignment=f"{PHYLO}/alignments/raw/{{marker}}.faa",
        raw_qc=f"{PHYLO}/alignments/raw/{{marker}}.json",
        code=f"{SCRIPTS}/infer_phylogeny.py",
        helpers=[f"{SCRIPTS}/busco_phylogeny.py", f"{SCRIPTS}/common.py"],
        binary=[PHY["trimal_command"]] if "/" in PHY["trimal_command"] else []
    output:
        alignment=f"{PHYLO}/alignments/{{marker}}.faa",
        qc=f"{PHYLO}/alignments/{{marker}}.json",
        columns=f"{PHYLO}/alignments/{{marker}}.columns.tsv"
    params:
        command=PHY["trimal_command"], mode=PHY["trimal_mode"],
        settings=json.dumps({k: PHY[k] for k in ["min_taxa", "min_protein_length"]}, sort_keys=True)
    threads: 1
    resources: mem_mb=PHY["trimming_mem_gb"] * 1000
    conda: "../envs/phylogeny.yaml"
    log: f"{LOG}/phylogeny/trim/{{marker}}.log"
    benchmark: f"{PHYLO}/benchmarks/trim.{{marker}}.tsv"
    shell:
        "{PYTHON:q} {input.code:q} trim --alignment {input.alignment:q} --raw-qc {input.raw_qc:q} "
        "--output {output.alignment:q} --qc {output.qc:q} --columns {output.columns:q} "
        "--command {params.command:q} --mode {params.mode:q} --settings {params.settings:q} > {log:q} 2>&1"


rule infer_busco_gene_tree:
    input:
        alignment=f"{PHYLO}/alignments/{{marker}}.faa",
        alignment_qc=f"{PHYLO}/alignments/{{marker}}.json",
        code=f"{SCRIPTS}/infer_phylogeny.py",
        helpers=[f"{SCRIPTS}/busco_phylogeny.py", f"{SCRIPTS}/common.py"]
    output:
        tree=f"{PHYLO}/gene_trees/{{marker}}.nwk",
        qc=f"{PHYLO}/gene_trees/{{marker}}.json"
    params: command=PHY["veryfasttree_command"], seed=PHY["seed"]
    threads: PHY["tree_threads"]
    resources: mem_mb=PHY["tree_mem_gb"] * 1000
    conda: "../envs/phylogeny.yaml"
    log: f"{LOG}/phylogeny/tree/{{marker}}.log"
    benchmark: f"{PHYLO}/benchmarks/tree.{{marker}}.tsv"
    shell:
        "{PYTHON:q} {input.code:q} gene_tree --alignment {input.alignment:q} --alignment-qc {input.alignment_qc:q} "
        "--output {output.tree:q} --qc {output.qc:q} --command {params.command:q} --threads {threads} "
        "--seed {params.seed} > {log:q} 2>&1"


rule merge_busco_gene_trees:
    input:
        manifest=f"{PHYLO_PLAN}/species.tsv",
        markers=f"{PHYLO_PLAN}/markers.tsv",
        trees=lambda wc: phylogeny_marker_files(wc, "nwk"),
        reports=lambda wc: phylogeny_marker_files(wc, "json"),
        code=f"{SCRIPTS}/infer_phylogeny.py",
        helpers=[f"{SCRIPTS}/busco_phylogeny.py", f"{SCRIPTS}/common.py"]
    output:
        trees=f"{PHYLO}/gene_trees.nwk",
        coverage=f"{PHYLO}/species_coverage.tsv",
        qc=f"{PHYLO}/gene_trees.json"
    params:
        tree_dir=f"{PHYLO}/gene_trees"
    conda: "../envs/phylogeny.yaml"
    resources: mem_mb=PHY["preparation_mem_gb"] * 1000
    log: f"{LOG}/phylogeny/merge.log"
    shell:
        "{PYTHON:q} {input.code:q} merge --manifest {input.manifest:q} --markers {input.markers:q} "
        "--tree-dir {params.tree_dir:q} --output {output.trees:q} --coverage {output.coverage:q} "
        "--qc {output.qc:q} > {log:q} 2>&1"


rule infer_busco_species_tree:
    input:
        trees=f"{PHYLO}/gene_trees.nwk",
        merge_qc=f"{PHYLO}/gene_trees.json",
        manifest=f"{PHYLO_PLAN}/species.tsv",
        code=f"{SCRIPTS}/infer_phylogeny.py",
        helpers=[f"{SCRIPTS}/busco_phylogeny.py", f"{SCRIPTS}/common.py", f"{SCRIPTS}/prepare_phylogeny_tools.py"],
        build_provenance=[str(Path(PHY["astral_command"]).parent.parent / "aster.json")] if "/" in PHY["astral_command"] and (Path(PHY["astral_command"]).parent.parent / "aster.json").is_file() else [],
        binary=[PHY["astral_command"]] if "/" in PHY["astral_command"] else []
    output: tree=f"{PHYLO}/species_tree.nwk", qc=f"{PHYLO}/species_tree.json"
    params: command=PHY["astral_command"], outgroup=PHY["outgroup"], seed=PHY["seed"]
    threads: PHY["astral_threads"]
    resources: mem_mb=PHY["astral_mem_gb"] * 1000
    conda: "../envs/phylogeny.yaml"
    log: f"{LOG}/phylogeny/astral.log"
    benchmark: f"{PHYLO}/benchmarks/astral.tsv"
    shell:
        "{PYTHON:q} {input.code:q} astral --trees {input.trees:q} --merge-qc {input.merge_qc:q} "
        "--manifest {input.manifest:q} --output {output.tree:q} --qc {output.qc:q} "
        "--command {params.command:q} --outgroup {params.outgroup:q} --threads {threads} "
        "--seed {params.seed} > {log:q} 2>&1"


rule prepare_timetree_calibrations:
    input:
        tree=f"{PHYLO}/species_tree.nwk",
        metadata=f"{META}/metadata_high_busco.tsv",
        taxonomy=config["taxonomy"]["database"],
        coverage=f"{PHYLO}/species_coverage.tsv",
        representatives=[PHY["dating"]["timetree"]["representatives"]] if PHY["dating"]["timetree"]["representatives"] else [],
        code=f"{SCRIPTS}/timetree_calibrations.py",
        helpers=[f"{SCRIPTS}/date_phylogeny.py", f"{SCRIPTS}/infer_phylogeny.py", f"{SCRIPTS}/busco_phylogeny.py", f"{SCRIPTS}/common.py"]
    output:
        calibrations=f"{PHYLO}/timetree/calibrations.tsv",
        candidates=f"{PHYLO}/timetree/candidates.tsv",
        details=f"{PHYLO}/timetree/candidates.json",
        provenance=f"{PHYLO}/timetree/provenance.json",
        representatives=f"{PHYLO}/timetree/representatives.txt",
        skeleton=f"{PHYLO}/timetree/representatives.nwk",
        taxa=f"{PHYLO}/timetree/taxa.tsv"
    params:
        outdir=f"{PHYLO}/timetree", cache=PHY["dating"]["timetree"]["cache_dir"],
        python=PHY["dating"]["timetree"]["python"],
        representatives_flag="--representatives" if PHY["dating"]["timetree"]["representatives"] else "",
        settings=json.dumps({k: PHY["dating"]["timetree"][k] for k in ["max_representatives", "max_queries",
            "min_studies", "offline", "request_delay_seconds"]}, sort_keys=True)
    threads: 1
    conda: "../envs/timetree.yaml"
    resources: mem_mb=PHY["dating"]["mem_gb"] * 1000
    log: f"{LOG}/phylogeny/timetree.log"
    benchmark: f"{PHYLO}/benchmarks/timetree.tsv"
    shell:
        "{params.python:q} {input.code:q} --tree {input.tree:q} --metadata {input.metadata:q} "
        "--taxonomy-db {input.taxonomy:q} --coverage {input.coverage:q} --outdir {params.outdir:q} "
        "--cache-dir {params.cache:q} --settings {params.settings:q} "
        "{params.representatives_flag} {input.representatives:q} > {log:q} 2>&1"


rule date_busco_species_tree:
    input:
        tree=f"{PHYLO}/species_tree.nwk",
        provenance=f"{PHYLO}/species_tree.json",
        calibrations=dating_calibrations,
        code=f"{SCRIPTS}/date_phylogeny.py",
        helpers=[f"{SCRIPTS}/infer_phylogeny.py", f"{SCRIPTS}/busco_phylogeny.py", f"{SCRIPTS}/common.py"],
        binary=[PHY["dating"]["command"]] if "/" in PHY["dating"]["command"] else []
    output:
        tree=f"{PHYLO}/dating/species_tree.dated.nwk",
        provenance=f"{PHYLO}/dating/provenance.json",
        ages=f"{PHYLO}/dating/node_ages.tsv",
        calibrations=f"{PHYLO}/dating/calibrations.resolved.tsv",
        raw=f"{PHYLO}/dating/treepl.dated.nwk",
        config=f"{PHYLO}/dating/treepl.config.txt",
        input_tree=f"{PHYLO}/dating/treepl.input.nwk",
        cv=f"{PHYLO}/dating/cross_validation.tsv",
        replicates=f"{PHYLO}/dating/optimization_replicates.tsv",
        adjustments=f"{PHYLO}/dating/branch_length_adjustments.tsv"
        # Keep treepl_runs as undeclared diagnostic artifacts: Snakemake must
        # not delete native logs and configurations when dating fails.
    params:
        outdir=f"{PHYLO}/dating", command=PHY["dating"]["command"],
        settings=json.dumps(PHY["dating"]["treepl"]), seed=PHY["seed"]
    conda: "../envs/phylogeny.yaml"
    threads: PHY["dating"]["threads"]
    resources: mem_mb=PHY["dating"]["mem_gb"] * 1000
    log: f"{LOG}/phylogeny/dating.log"
    benchmark: f"{PHYLO}/benchmarks/dating.tsv"
    shell:
        "{PYTHON:q} {input.code:q} --tree {input.tree:q} --provenance {input.provenance:q} "
        "--calibrations {input.calibrations:q} --outdir {params.outdir:q} --command {params.command:q} "
        "--settings {params.settings:q} --threads {threads} --seed {params.seed} > {log:q} 2>&1"
