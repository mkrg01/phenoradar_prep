def phylogeny_samples(wc):
    if wc.phylo_branch == "contrast/phylogeny":
        return checkpoints.select_contrast_representatives.get().output.samples
    if wc.phylo_branch == "phylogeny_phenotyped":
        return checkpoints.select_phenotyped_species.get().output.samples
    return checkpoints.select_metadata.get().output.samples


def phylogeny_root_guide(wc):
    if PHY["outgroup"] != "auto":
        return []
    if wc.phylo_branch == "contrast/phylogeny":
        return [checkpoints.select_contrast_representatives.get().output.tree]
    if wc.phylo_branch == "phylogeny_phenotyped":
        return [f"{PHENOTYPED}/rooting/ncbi_tree.nwk"]
    return [f"{ROOTING}/ncbi_tree.nwk"]


def phylogeny_input_rows(wc):
    with open(phylogeny_samples(wc)) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def phylogeny_tables(wc):
    if not PHY["busco_full_dir"]:
        raise WorkflowError("set phylogeny.busco_full_dir to per-species BUSCO full tables; see docs/phylogeny.md")
    return sorted({str(Path(PHY["busco_full_dir"]) / (r["species"] + PHY["busco_full_suffix"])) for r in phylogeny_input_rows(wc)})


def phylogeny_plan_rows(wc, table):
    output = checkpoints.plan_phylogeny.get(phylo_branch=wc.phylo_branch).output
    with open(getattr(output, table)) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def phylogeny_species_row(wc):
    rows = [r for r in phylogeny_plan_rows(wc, "species") if r["species"] == wc.species]
    if len(rows) != 1:
        raise WorkflowError(f"species absent from phylogeny plan: {wc.species}")
    return rows[0]


def phylogeny_marker_files(wc, suffix):
    return [f'{OUT}/{wc.phylo_branch}/gene_trees/{r["marker"]}.{suffix}' for r in phylogeny_plan_rows(wc, "markers")]


def dating_calibrations(wc):
    if PHY["dating"]["calibration_source"] == "timetree":
        return f"{OUT}/{wc.phylo_branch}/timetree/calibrations.tsv"
    if not PHY["dating"]["calibrations"]:
        raise WorkflowError("absolute dating requires phylogeny.dating.calibrations; no arbitrary root age is used")
    return PHY["dating"]["calibrations"]


checkpoint select_phenotyped_species:
    input:
        samples=f"{META}/samples.tsv",
        traits=config["inputs"]["species_trait"],
        code=f"{SCRIPTS}/species_traits.py",
        common=f"{SCRIPTS}/common.py"
    output:
        samples=f"{PHENOTYPED}/selection/samples.tsv",
        qc=f"{PHENOTYPED}/selection/selection.json"
    params:
        outdir=f"{PHENOTYPED}/selection", trait=PHY["trait"], min_taxa=PHY["min_taxa"]
    conda: "../envs/analysis.yaml"
    resources: mem_mb=4000
    log: f"{LOG}/phylogeny_phenotyped/selection.log"
    shell:
        "{PYTHON:q} {input.code:q} --samples {input.samples:q} --traits {input.traits:q} "
        "--trait {params.trait:q} --min-taxa {params.min_taxa} --outdir {params.outdir:q} > {log:q} 2>&1"


checkpoint plan_phylogeny:
    input:
        samples=phylogeny_samples,
        outgroup=f"{PHYLO_RUN}/rooting/outgroup.txt",
        tables=phylogeny_tables,
        code=f"{SCRIPTS}/busco_phylogeny.py",
        common=f"{SCRIPTS}/common.py"
    output:
        species=f"{PHYLO_RUN}/plan/species.tsv",
        markers=f"{PHYLO_RUN}/plan/markers.tsv",
        stats=f"{PHYLO_RUN}/plan/marker_stats.tsv",
        provenance=f"{PHYLO_RUN}/plan/provenance.json"
    params:
        outdir=f"{PHYLO_RUN}/plan",
        settings=json.dumps({k: PHY[k] for k in ["outgroup", "busco_full_dir", "busco_full_suffix", "lineage",
            "sequence_dir", "sequence_suffix", "min_taxa", "max_markers"]}, sort_keys=True)
    log: f"{LOG}/{{phylo_branch}}/plan.log"
    conda: "../envs/phylogeny.yaml"
    resources: mem_mb=PHY["preparation_mem_gb"] * 1000
    shell:
        "{PYTHON:q} {input.code:q} plan --samples {input.samples:q} "
        "--outgroup-file {input.outgroup:q} --outdir {params.outdir:q} --settings {params.settings:q} > {log:q} 2>&1"


rule extract_busco_proteins:
    input:
        markers=f"{PHYLO_RUN}/plan/markers.tsv",
        table=lambda wc: phylogeny_species_row(wc)["busco_table"],
        sequences=lambda wc: phylogeny_species_row(wc)["sequences"],
        code=f"{SCRIPTS}/busco_phylogeny.py",
        common=f"{SCRIPTS}/common.py"
    output:
        proteins=f"{PHYLO_RUN}/species/{{species}}.faa",
        qc=f"{PHYLO_RUN}/species/{{species}}.json"
    params:
        settings=json.dumps({k: PHY[k] for k in ["lineage", "sequence_mode", "translation_table",
            "min_protein_length", "max_unknown_fraction"]}, sort_keys=True)
    conda: "../envs/phylogeny.yaml"
    resources: mem_mb=PHY["preparation_mem_gb"] * 1000
    log: f"{LOG}/{{phylo_branch}}/extract/{{species}}.log"
    shell:
        "{PYTHON:q} {input.code:q} extract --species {wildcards.species:q} --table {input.table:q} "
        "--sequences {input.sequences:q} --markers {input.markers:q} --output {output.proteins:q} "
        "--qc {output.qc:q} --settings {params.settings:q} > {log:q} 2>&1"


rule collect_busco_markers:
    input:
        manifest=f"{PHYLO_RUN}/plan/species.tsv",
        markers=f"{PHYLO_RUN}/plan/markers.tsv",
        proteins=lambda wc: [f'{OUT}/{wc.phylo_branch}/species/{r["species"]}.faa' for r in phylogeny_plan_rows(wc, "species")],
        reports=lambda wc: [f'{OUT}/{wc.phylo_branch}/species/{r["species"]}.json' for r in phylogeny_plan_rows(wc, "species")],
        code=f"{SCRIPTS}/busco_phylogeny.py",
        common=f"{SCRIPTS}/common.py"
    output: fasta=directory(f"{PHYLO_RUN}/markers")
    params: species_dir=f"{PHYLO_RUN}/species"
    conda: "../envs/phylogeny.yaml"
    resources: mem_mb=PHY["preparation_mem_gb"] * 1000
    log: f"{LOG}/{{phylo_branch}}/collect.log"
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
        alignment=f"{PHYLO_RUN}/alignments/raw/{{marker}}.faa",
        qc=f"{PHYLO_RUN}/alignments/raw/{{marker}}.json"
    params:
        fasta=lambda wc: f"{OUT}/{wc.phylo_branch}/markers/{wc.marker}.faa",
        command="famsa",
        settings=json.dumps({"min_taxa": PHY["min_taxa"]}, sort_keys=True)
    threads: PHY["align_threads"]
    resources: mem_mb=PHY["alignment_mem_gb"] * 1000
    conda: "../envs/phylogeny.yaml"
    log: f"{LOG}/{{phylo_branch}}/align/{{marker}}.log"
    benchmark: f"{PHYLO_RUN}/benchmarks/align.{{marker}}.tsv"
    shell:
        "{PYTHON:q} {input.code:q} align --fasta {params.fasta:q} --output {output.alignment:q} "
        "--qc {output.qc:q} --command {params.command:q} --threads {threads} "
        "--settings {params.settings:q} > {log:q} 2>&1"


rule trim_busco_marker:
    wildcard_constraints: marker="[A-Za-z0-9][A-Za-z0-9_.-]*"
    input:
        alignment=f"{PHYLO_RUN}/alignments/raw/{{marker}}.faa",
        raw_qc=f"{PHYLO_RUN}/alignments/raw/{{marker}}.json",
        code=f"{SCRIPTS}/infer_phylogeny.py",
        helpers=[f"{SCRIPTS}/busco_phylogeny.py", f"{SCRIPTS}/common.py"]
    output:
        alignment=f"{PHYLO_RUN}/alignments/{{marker}}.faa",
        qc=f"{PHYLO_RUN}/alignments/{{marker}}.json",
        columns=f"{PHYLO_RUN}/alignments/{{marker}}.columns.tsv"
    params:
        command="trimal", mode=PHY["trimal_mode"],
        settings=json.dumps({k: PHY[k] for k in ["min_taxa", "min_protein_length"]}, sort_keys=True)
    threads: 1
    resources: mem_mb=PHY["trimming_mem_gb"] * 1000
    conda: "../envs/phylogeny.yaml"
    log: f"{LOG}/{{phylo_branch}}/trim/{{marker}}.log"
    benchmark: f"{PHYLO_RUN}/benchmarks/trim.{{marker}}.tsv"
    shell:
        "{PYTHON:q} {input.code:q} trim --alignment {input.alignment:q} --raw-qc {input.raw_qc:q} "
        "--output {output.alignment:q} --qc {output.qc:q} --columns {output.columns:q} "
        "--command {params.command:q} --mode {params.mode:q} --settings {params.settings:q} > {log:q} 2>&1"


rule infer_busco_gene_tree:
    input:
        alignment=f"{PHYLO_RUN}/alignments/{{marker}}.faa",
        alignment_qc=f"{PHYLO_RUN}/alignments/{{marker}}.json",
        code=f"{SCRIPTS}/infer_phylogeny.py",
        helpers=[f"{SCRIPTS}/busco_phylogeny.py", f"{SCRIPTS}/common.py"]
    output:
        tree=f"{PHYLO_RUN}/gene_trees/{{marker}}.nwk",
        qc=f"{PHYLO_RUN}/gene_trees/{{marker}}.json"
    params: command="VeryFastTree", seed=PHY["seed"]
    threads: PHY["tree_threads"]
    resources: mem_mb=PHY["tree_mem_gb"] * 1000
    conda: "../envs/phylogeny.yaml"
    log: f"{LOG}/{{phylo_branch}}/tree/{{marker}}.log"
    benchmark: f"{PHYLO_RUN}/benchmarks/tree.{{marker}}.tsv"
    shell:
        "{PYTHON:q} {input.code:q} gene_tree --alignment {input.alignment:q} --alignment-qc {input.alignment_qc:q} "
        "--output {output.tree:q} --qc {output.qc:q} --command {params.command:q} --threads {threads} "
        "--seed {params.seed} > {log:q} 2>&1"


rule merge_busco_gene_trees:
    input:
        manifest=f"{PHYLO_RUN}/plan/species.tsv",
        markers=f"{PHYLO_RUN}/plan/markers.tsv",
        trees=lambda wc: phylogeny_marker_files(wc, "nwk"),
        reports=lambda wc: phylogeny_marker_files(wc, "json"),
        code=f"{SCRIPTS}/infer_phylogeny.py",
        helpers=[f"{SCRIPTS}/busco_phylogeny.py", f"{SCRIPTS}/common.py"]
    output:
        trees=f"{PHYLO_RUN}/gene_trees.nwk",
        coverage=f"{PHYLO_RUN}/species_coverage.tsv",
        qc=f"{PHYLO_RUN}/gene_trees.json"
    params:
        tree_dir=f"{PHYLO_RUN}/gene_trees"
    conda: "../envs/phylogeny.yaml"
    resources: mem_mb=PHY["preparation_mem_gb"] * 1000
    log: f"{LOG}/{{phylo_branch}}/merge.log"
    shell:
        "{PYTHON:q} {input.code:q} merge --manifest {input.manifest:q} --markers {input.markers:q} "
        "--tree-dir {params.tree_dir:q} --output {output.trees:q} --coverage {output.coverage:q} "
        "--qc {output.qc:q} > {log:q} 2>&1"


rule infer_busco_species_tree:
    input:
        trees=f"{PHYLO_RUN}/gene_trees.nwk",
        merge_qc=f"{PHYLO_RUN}/gene_trees.json",
        manifest=f"{PHYLO_RUN}/plan/species.tsv",
        code=f"{SCRIPTS}/infer_phylogeny.py",
        helpers=[f"{SCRIPTS}/busco_phylogeny.py", f"{SCRIPTS}/common.py", f"{SCRIPTS}/prepare_phylogeny_tools.py"],
        build_provenance=str(Path(ASTRAL).parent.parent / "aster.json"),
        binary=ASTRAL,
        outgroup=f"{PHYLO_RUN}/rooting/outgroup.txt"
    output: tree=f"{PHYLO_RUN}/species_tree.nwk", qc=f"{PHYLO_RUN}/species_tree.json"
    params: command=ASTRAL, seed=PHY["seed"]
    threads: PHY["astral_threads"]
    resources: mem_mb=PHY["astral_mem_gb"] * 1000
    conda: "../envs/phylogeny.yaml"
    log: f"{LOG}/{{phylo_branch}}/astral.log"
    benchmark: f"{PHYLO_RUN}/benchmarks/astral.tsv"
    shell:
        "{PYTHON:q} {input.code:q} astral --trees {input.trees:q} --merge-qc {input.merge_qc:q} "
        "--manifest {input.manifest:q} --output {output.tree:q} --qc {output.qc:q} "
        "--command {params.command:q} --outgroup-file {input.outgroup:q} --threads {threads} "
        "--seed {params.seed} > {log:q} 2>&1"


rule prepare_timetree_calibrations:
    wildcard_constraints: phylo_branch="phylogeny|phylogeny_phenotyped"
    input:
        tree=f"{PHYLO_RUN}/species_tree.nwk",
        metadata=f"{META}/metadata_high_busco.tsv",
        taxonomy=TAXONOMY_DB,
        coverage=f"{PHYLO_RUN}/species_coverage.tsv",
        representatives=[PHY["dating"]["timetree"]["representatives"]] if PHY["dating"]["timetree"]["representatives"] else [],
        code=f"{SCRIPTS}/timetree_calibrations.py",
        helpers=[f"{SCRIPTS}/date_phylogeny.py", f"{SCRIPTS}/infer_phylogeny.py", f"{SCRIPTS}/busco_phylogeny.py", f"{SCRIPTS}/common.py"]
    output:
        calibrations=f"{PHYLO_RUN}/timetree/calibrations.tsv",
        candidates=f"{PHYLO_RUN}/timetree/candidates.tsv",
        details=f"{PHYLO_RUN}/timetree/candidates.json",
        provenance=f"{PHYLO_RUN}/timetree/provenance.json",
        representatives=f"{PHYLO_RUN}/timetree/representatives.txt",
        skeleton=f"{PHYLO_RUN}/timetree/representatives.nwk",
        taxa=f"{PHYLO_RUN}/timetree/taxa.tsv"
    params:
        outdir=f"{PHYLO_RUN}/timetree", cache=TIMETREE_CACHE,
        representatives_flag="--representatives" if PHY["dating"]["timetree"]["representatives"] else "",
        settings=json.dumps(dict({k: PHY["dating"]["timetree"][k] for k in ["max_representatives", "max_queries",
            "min_studies", "offline"]}, request_delay_seconds=1.0), sort_keys=True)
    threads: 1
    conda: "../envs/timetree.yaml"
    resources: mem_mb=PHY["dating"]["mem_gb"] * 1000
    log: f"{LOG}/{{phylo_branch}}/timetree.log"
    benchmark: f"{PHYLO_RUN}/benchmarks/timetree.tsv"
    shell:
        "{PYTHON:q} {input.code:q} --tree {input.tree:q} --metadata {input.metadata:q} "
        "--taxonomy-db {input.taxonomy:q} --coverage {input.coverage:q} --outdir {params.outdir:q} "
        "--cache-dir {params.cache:q} --settings {params.settings:q} "
        "{params.representatives_flag} {input.representatives:q} > {log:q} 2>&1"


rule date_busco_species_tree:
    wildcard_constraints: phylo_branch="phylogeny|phylogeny_phenotyped"
    input:
        tree=f"{PHYLO_RUN}/species_tree.nwk",
        provenance=f"{PHYLO_RUN}/species_tree.json",
        calibrations=dating_calibrations,
        code=f"{SCRIPTS}/date_phylogeny.py",
        helpers=[f"{SCRIPTS}/infer_phylogeny.py", f"{SCRIPTS}/busco_phylogeny.py", f"{SCRIPTS}/common.py"]
    output:
        tree=f"{PHYLO_RUN}/dating/species_tree.dated.nwk",
        provenance=f"{PHYLO_RUN}/dating/provenance.json",
        ages=f"{PHYLO_RUN}/dating/node_ages.tsv",
        calibrations=f"{PHYLO_RUN}/dating/calibrations.resolved.tsv",
        raw=f"{PHYLO_RUN}/dating/lsd2.dated.date.nexus",
        fitted=f"{PHYLO_RUN}/dating/lsd2.fitted.nwk",
        report=f"{PHYLO_RUN}/dating/lsd2.report.txt",
        dates=f"{PHYLO_RUN}/dating/lsd2.dates.txt",
        command=f"{PHYLO_RUN}/dating/lsd2.command.json",
        input_tree=f"{PHYLO_RUN}/dating/lsd2.input.nwk",
        adjustments=f"{PHYLO_RUN}/dating/rounding_adjustments.tsv"
        # Keep lsd2_runs undeclared so native diagnostics survive a failed job.
    params:
        outdir=f"{PHYLO_RUN}/dating", command="lsd2",
        settings=json.dumps(PHY["dating"]["lsd2"])
    conda: "../envs/dating.yaml"
    threads: 1
    resources: mem_mb=PHY["dating"]["mem_gb"] * 1000
    log: f"{LOG}/{{phylo_branch}}/dating.log"
    benchmark: f"{PHYLO_RUN}/benchmarks/dating.tsv"
    shell:
        "{PYTHON:q} {input.code:q} --tree {input.tree:q} --provenance {input.provenance:q} "
        "--calibrations {input.calibrations:q} --outdir {params.outdir:q} --command {params.command:q} "
        "--settings {params.settings:q} --threads {threads} > {log:q} 2>&1"
