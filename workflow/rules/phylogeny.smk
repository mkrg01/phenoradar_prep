from busco_phylogeny import busco_full_path


def phylogeny_samples(wc):
    report_phylogeny_plan()
    if wc.phylo_branch == REPRESENTATIVES_REL:
        return checkpoints.select_contrast_representatives.get().output.samples
    if wc.phylo_branch == PHYLO_BRANCHES["phenotyped"]:
        return checkpoints.select_phenotyped_species.get().output.samples
    return checkpoints.select_metadata.get().output.samples


def phylogeny_root_guide(wc):
    if PHY["outgroup"] != "auto":
        return []
    if wc.phylo_branch == REPRESENTATIVES_REL:
        return [checkpoints.select_contrast_representatives.get().output.tree]
    if wc.phylo_branch == PHYLO_BRANCHES["phenotyped"]:
        return [f"{PHENOTYPED}/rooting/ncbi_tree.nwk"]
    return [f"{ROOTING}/ncbi_tree.nwk"]


def phylogeny_input_rows(wc):
    with open(phylogeny_samples(wc)) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def phylogeny_tables(wc):
    return [str(busco_full_path(BUSCO_FULL, species, PHY["lineage"]))
            for species in sorted({r["species"] for r in phylogeny_input_rows(wc)})]


REPORTED_PHYLOGENY_PLANS = set()


def phylogeny_plan_rows(wc, table):
    output = checkpoints.plan_phylogeny.get(phylo_branch=wc.phylo_branch).output
    if wc.phylo_branch not in REPORTED_PHYLOGENY_PLANS:
        with open(output.species) as handle:
            count = sum(1 for _ in csv.DictReader(handle, delimiter="\t"))
        logger.info(f"Phylogeny inference species: {wc.phylo_branch} = {count} (trait={config['trait']})")
        REPORTED_PHYLOGENY_PLANS.add(wc.phylo_branch)
    with open(getattr(output, table)) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def phylogeny_species_row(wc):
    rows = [r for r in phylogeny_plan_rows(wc, "species") if r["species"] == wc.species]
    if len(rows) != 1:
        raise WorkflowError(f"species absent from phylogeny plan: {wc.species}")
    return rows[0]


def phylogeny_marker_files(wc, suffix):
    return [f'{OUT}/{wc.phylo_branch}/gene_trees/{r["marker"]}.{suffix}' for r in phylogeny_plan_rows(wc, "markers")]


checkpoint select_phenotyped_species:
    input:
        samples=f"{META}/samples.tsv",
        traits=rules.normalize_traits.output.traits,
        code=f"{SCRIPTS}/species_traits.py",
        common=[f"{SCRIPTS}/common.py", f"{SCRIPTS}/sample_identity.py", f"{SCRIPTS}/phylogeny_outgroup.py"]
    output:
        samples=f"{PHENOTYPED}/selection/samples.tsv",
        qc=f"{PHENOTYPED}/selection/selection.json"
    params:
        outdir=f"{PHENOTYPED}/selection", trait=PHY["trait"]
    conda: "../envs/analysis.yaml"
    resources: mem_mb=4000
    log: f"{LOG}/{PHYLO_BRANCHES['phenotyped']}/selection.log"
    shell:
        "{PYTHON:q} {input.code:q} --samples {input.samples:q} --traits {input.traits:q} "
        "--trait {params.trait:q} --outdir {params.outdir:q} > {log:q} 2>&1"


checkpoint plan_phylogeny:
    input:
        samples=phylogeny_samples,
        anchor=f"{PHYLO_RUN}/rooting/anchor.txt",
        tables=phylogeny_tables,
        code=f"{SCRIPTS}/busco_phylogeny.py",
        common=[f"{SCRIPTS}/common.py", f"{SCRIPTS}/sample_identity.py", f"{SCRIPTS}/phylogeny_outgroup.py"]
    output:
        species=f"{PHYLO_RUN}/plan/species.tsv",
        markers=f"{PHYLO_RUN}/plan/markers.tsv",
        stats=f"{PHYLO_RUN}/plan/marker_stats.tsv",
        provenance=f"{PHYLO_RUN}/plan/provenance.json"
    params:
        outdir=f"{PHYLO_RUN}/plan",
        settings=json.dumps({k: PHY[k] for k in ["outgroup", "busco_full_dir", "lineage",
            "min_taxa", "max_markers"]}, sort_keys=True)
    log: f"{LOG}/{{phylo_branch}}/plan.log"
    conda: "../envs/phylogeny.yaml"
    resources: mem_mb=4000
    shell:
        "{PYTHON:q} {input.code:q} plan --samples {input.samples:q} "
        "--outgroup-file {input.anchor:q} --outdir {params.outdir:q} --settings {params.settings:q} > {log:q} 2>&1"


rule extract_busco_proteins:
    input:
        markers=f"{PHYLO_RUN}/plan/markers.tsv",
        table=lambda wc: phylogeny_species_row(wc)["busco_table"],
        sequences=lambda wc: phylogeny_species_row(wc)["sequences"],
        code=f"{SCRIPTS}/busco_phylogeny.py",
        common=[f"{SCRIPTS}/common.py", f"{SCRIPTS}/sample_identity.py", f"{SCRIPTS}/phylogeny_outgroup.py"]
    output:
        proteins=f"{PHYLO_RUN}/species/{{species}}.faa",
        qc=f"{PHYLO_RUN}/species/{{species}}.json"
    params:
        settings=json.dumps({k: PHY[k] for k in ["lineage", "translation_table",
            "min_protein_length", "max_unknown_fraction"]}, sort_keys=True)
    conda: "../envs/phylogeny.yaml"
    resources: mem_mb=4000
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
        common=[f"{SCRIPTS}/common.py", f"{SCRIPTS}/sample_identity.py", f"{SCRIPTS}/phylogeny_outgroup.py"]
    output: fasta=directory(f"{PHYLO_RUN}/markers")
    params: species_dir=f"{PHYLO_RUN}/species"
    conda: "../envs/phylogeny.yaml"
    resources: mem_mb=4000
    log: f"{LOG}/{{phylo_branch}}/collect.log"
    shell:
        "{PYTHON:q} {input.code:q} collect --manifest {input.manifest:q} --markers {input.markers:q} "
        "--species-dir {params.species_dir:q} --outdir {output.fasta:q} > {log:q} 2>&1"


def busco_marker_mem_mb(wildcards):
    return 4000 if wildcards.phylo_branch == REPRESENTATIVES_REL else 8000


# Keep each marker chain in one cluster allocation; completed stages remain reusable.
rule align_busco_marker:
    group: "busco_marker"
    wildcard_constraints: marker="[A-Za-z0-9][A-Za-z0-9_.-]*"
    input:
        markers=rules.collect_busco_markers.output.fasta,
        sample_manifest=f"{PHYLO_RUN}/plan/species.tsv",
        code=f"{SCRIPTS}/infer_phylogeny.py",
        helpers=[f"{SCRIPTS}/busco_phylogeny.py", f"{SCRIPTS}/common.py", f"{SCRIPTS}/sample_identity.py", f"{SCRIPTS}/phylogeny_outgroup.py"]
    output:
        alignment=f"{PHYLO_RUN}/alignments/raw/{{marker}}.faa",
        qc=f"{PHYLO_RUN}/alignments/raw/{{marker}}.json"
    params:
        fasta=lambda wc: f"{OUT}/{wc.phylo_branch}/markers/{wc.marker}.faa",
        command="famsa",
        settings=lambda wc: json.dumps({"min_taxa": PHY["min_taxa"], "sample_manifest": f"{OUT}/{wc.phylo_branch}/plan/species.tsv"}, sort_keys=True)
    threads: 1
    resources: mem_mb=16000
    conda: "../envs/phylogeny.yaml"
    log: f"{LOG}/{{phylo_branch}}/align/{{marker}}.log"
    benchmark: f"{LOG}/{{phylo_branch}}/benchmarks/align.{{marker}}.tsv"
    shell:
        "{PYTHON:q} {input.code:q} align --fasta {params.fasta:q} --output {output.alignment:q} "
        "--qc {output.qc:q} --command {params.command:q} --threads {threads} "
        "--settings {params.settings:q} > {log:q} 2>&1"


rule trim_busco_marker:
    group: "busco_marker"
    wildcard_constraints: marker="[A-Za-z0-9][A-Za-z0-9_.-]*"
    input:
        alignment=f"{PHYLO_RUN}/alignments/raw/{{marker}}.faa",
        raw_qc=f"{PHYLO_RUN}/alignments/raw/{{marker}}.json",
        sample_manifest=f"{PHYLO_RUN}/plan/species.tsv",
        code=f"{SCRIPTS}/infer_phylogeny.py",
        helpers=[f"{SCRIPTS}/busco_phylogeny.py", f"{SCRIPTS}/common.py", f"{SCRIPTS}/sample_identity.py", f"{SCRIPTS}/phylogeny_outgroup.py"],
        trimming=f"{SCRIPTS}/alignment_trimming.py"
    output:
        alignment=f"{PHYLO_RUN}/alignments/{{marker}}.faa",
        qc=f"{PHYLO_RUN}/alignments/{{marker}}.json",
        columns=f"{PHYLO_RUN}/alignments/{{marker}}.columns.tsv"
    params:
        command="trimal", mode=PHY["trimal_mode"],
        settings=lambda wc: json.dumps({**{k: PHY[k] for k in ["min_taxa", "min_protein_length"]}, "sample_manifest": f"{OUT}/{wc.phylo_branch}/plan/species.tsv"}, sort_keys=True)
    threads: 1
    resources: mem_mb=4000
    conda: "../envs/phylogeny.yaml"
    log: f"{LOG}/{{phylo_branch}}/trim/{{marker}}.log"
    benchmark: f"{LOG}/{{phylo_branch}}/benchmarks/trim.{{marker}}.tsv"
    shell:
        "{PYTHON:q} {input.code:q} trim --alignment {input.alignment:q} --raw-qc {input.raw_qc:q} "
        "--output {output.alignment:q} --qc {output.qc:q} --columns {output.columns:q} "
        "--command {params.command:q} --mode {params.mode:q} --settings {params.settings:q} > {log:q} 2>&1"


rule infer_busco_gene_tree:
    group: "busco_marker"
    input:
        alignment=f"{PHYLO_RUN}/alignments/{{marker}}.faa",
        alignment_qc=f"{PHYLO_RUN}/alignments/{{marker}}.json",
        code=f"{SCRIPTS}/infer_phylogeny.py",
        helpers=[f"{SCRIPTS}/busco_phylogeny.py", f"{SCRIPTS}/common.py", f"{SCRIPTS}/sample_identity.py", f"{SCRIPTS}/phylogeny_outgroup.py"]
    output:
        tree=f"{PHYLO_RUN}/gene_trees/{{marker}}.nwk",
        qc=f"{PHYLO_RUN}/gene_trees/{{marker}}.json"
    params: command="VeryFastTree", seed=RANDOM_SEED
    threads: 1
    resources: mem_mb=busco_marker_mem_mb
    conda: "../envs/phylogeny.yaml"
    log: f"{LOG}/{{phylo_branch}}/tree/{{marker}}.log"
    benchmark: f"{LOG}/{{phylo_branch}}/benchmarks/tree.{{marker}}.tsv"
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
        helpers=[f"{SCRIPTS}/busco_phylogeny.py", f"{SCRIPTS}/common.py", f"{SCRIPTS}/sample_identity.py", f"{SCRIPTS}/phylogeny_outgroup.py"]
    output:
        trees=f"{PHYLO_RUN}/gene_trees.nwk",
        coverage=f"{PHYLO_RUN}/species_coverage.tsv",
        qc=f"{PHYLO_RUN}/gene_trees.json"
    params:
        tree_dir=f"{PHYLO_RUN}/gene_trees"
    conda: "../envs/phylogeny.yaml"
    resources: mem_mb=4000
    log: f"{LOG}/{{phylo_branch}}/merge.log"
    shell:
        "{PYTHON:q} {input.code:q} merge --manifest {input.manifest:q} --markers {input.markers:q} "
        "--tree-dir {params.tree_dir:q} --output {output.trees:q} --coverage {output.coverage:q} "
        "--qc {output.qc:q} > {log:q} 2>&1"


# Native Conda builds the pinned official source; containers export the verified
# binary built into the image. Existing user-prepared tools remain external inputs.
if not Path(ASTRAL).is_file() or not Path(ASTRAL).parent.parent.joinpath("aster.json").is_file():
    rule prepare_astral:
        input:
            code=f"{SCRIPTS}/prepare_phylogeny_tools.py",
            common=[f"{SCRIPTS}/common.py", f"{SCRIPTS}/sample_identity.py", f"{SCRIPTS}/phylogeny_outgroup.py"]
        output:
            binary=ASTRAL,
            provenance=str(Path(ASTRAL).parent.parent / "aster.json")
        params: destination=str(Path(ASTRAL).parent.parent)
        conda: "../envs/dating.yaml"
        threads: 1
        resources: mem_mb=8000
        log: f"{LOG}/astral_prepare.log"
        shell:
            "{PYTHON:q} {input.code:q} --destination {params.destination:q} > {log:q} 2>&1"


rule infer_busco_species_tree:
    input:
        trees=f"{PHYLO_RUN}/gene_trees.nwk",
        merge_qc=f"{PHYLO_RUN}/gene_trees.json",
        manifest=f"{PHYLO_RUN}/plan/species.tsv",
        code=f"{SCRIPTS}/infer_phylogeny.py",
        helpers=[f"{SCRIPTS}/busco_phylogeny.py", f"{SCRIPTS}/common.py", f"{SCRIPTS}/prepare_phylogeny_tools.py", f"{SCRIPTS}/phylogeny_outgroup.py"],
        build_provenance=str(Path(ASTRAL).parent.parent / "aster.json"),
        binary=ASTRAL,
        anchor=f"{PHYLO_RUN}/rooting/anchor.txt"
    output:
        tree=f"{PHYLO_RUN}/rooting/astral.nwk",
        qc=f"{PHYLO_RUN}/rooting/astral.json"
    params: command=ASTRAL, seed=RANDOM_SEED
    threads: 32
    resources: mem_mb=128000
    conda: "../envs/phylogeny.yaml"
    log: f"{LOG}/{{phylo_branch}}/astral.log"
    benchmark: f"{LOG}/{{phylo_branch}}/benchmarks/astral.tsv"
    shell:
        "{PYTHON:q} {input.code:q} astral --trees {input.trees:q} --merge-qc {input.merge_qc:q} "
        "--manifest {input.manifest:q} --output {output.tree:q} --qc {output.qc:q} "
        "--command {params.command:q} --outgroup-file {input.anchor:q} --threads {threads} "
        "--seed {params.seed} > {log:q} 2>&1"


rule root_busco_species_tree:
    input:
        samples=phylogeny_samples,
        taxonomy=TAXONOMY_DB,
        anchor_qc=f"{PHYLO_RUN}/rooting/anchor.json",
        tree=f"{PHYLO_RUN}/rooting/astral.nwk",
        tree_qc=f"{PHYLO_RUN}/rooting/astral.json",
        code=f"{SCRIPTS}/phylogeny_root.py",
        helpers=[f"{SCRIPTS}/{name}.py" for name in
                 ["common", "phylogeny_outgroup", "sample_identity", "infer_phylogeny", "busco_phylogeny", "nwkit_runtime"]]
    output:
        tree=f"{PHYLO_RUN}/species_tree.nwk",
        qc=f"{PHYLO_RUN}/species_tree.json",
        outgroup=f"{PHYLO_RUN}/rooting/outgroup.txt",
        rooting_qc=f"{PHYLO_RUN}/rooting/outgroup.json"
    conda: "../envs/dating.yaml"
    resources: mem_mb=4000
    log: f"{LOG}/{{phylo_branch}}/rooting/outgroup.log"
    shell:
        "{PYTHON:q} {input.code:q} root_species_tree --samples {input.samples:q} --taxonomy-db {input.taxonomy:q} "
        "--tree {input.tree:q} --tree-qc {input.tree_qc:q} --anchor-qc {input.anchor_qc:q} "
        "--output {output.tree:q} --qc {output.qc:q} --outgroup-file {output.outgroup:q} "
        "--rooting-qc {output.rooting_qc:q} > {log:q} 2>&1"


def dating_samples(wc):
    if wc.phylo_branch == PHYLO_BRANCHES["phenotyped"]:
        return f"{PHENOTYPED}/selection/samples.tsv"
    return f"{META}/samples.tsv"


rule prepare_angiocal_calibrations:
    wildcard_constraints: phylo_branch=DATING_BRANCH_PATTERN
    input:
        tree=f"{PHYLO_RUN}/species_tree.nwk",
        samples=dating_samples,
        taxonomy=TAXONOMY_DB,
        mapping=[FOSSIL_MAP] if Path(FOSSIL_MAP).is_file() else [],
        code=f"{SCRIPTS}/angiocal_calibrations.py",
        helpers=[f"{SCRIPTS}/{name}.py" for name in
                 ["infer_phylogeny", "taxonomy_check", "busco_phylogeny",
                  "common", "sample_identity", "phylogeny_outgroup", "nwkit_runtime"]]
    output:
        calibrations=f"{PHYLO_RUN}/angiocal/calibrations.tsv",
        fossils=f"{PHYLO_RUN}/angiocal/fossils.tsv",
        provenance=f"{PHYLO_RUN}/angiocal/provenance.json"
    params: outdir=f"{PHYLO_RUN}/angiocal", map_flag="--mapping" if Path(FOSSIL_MAP).is_file() else ""
    threads: 1
    conda: "../envs/dating.yaml"
    resources: mem_mb=4000
    log: f"{LOG}/{{phylo_branch}}/angiocal.log"
    benchmark: f"{LOG}/{{phylo_branch}}/benchmarks/angiocal.tsv"
    shell:
        "{PYTHON:q} {input.code:q} --tree {input.tree:q} --samples {input.samples:q} "
        "--taxonomy {input.taxonomy:q} --outdir {params.outdir:q} "
        "{params.map_flag} {input.mapping:q} > {log:q} 2>&1"


rule prepare_timetree_calibrations:
    wildcard_constraints: phylo_branch=DATING_BRANCH_PATTERN
    input:
        tree=f"{PHYLO_RUN}/species_tree.nwk",
        samples=dating_samples,
        taxonomy=TAXONOMY_DB,
        code=f"{SCRIPTS}/timetree_calibrations.py",
        helpers=[f"{SCRIPTS}/{name}.py" for name in
                 ["infer_phylogeny", "taxonomy_check", "busco_phylogeny",
                  "common", "sample_identity", "phylogeny_outgroup", "nwkit_runtime"]]
    output:
        calibrations=f"{PHYLO_RUN}/timetree/calibrations.tsv",
        provenance=f"{PHYLO_RUN}/timetree/provenance.json"
    params: outdir=f"{PHYLO_RUN}/timetree"
    threads: 1
    conda: "../envs/dating.yaml"
    resources: mem_mb=4000
    log: f"{LOG}/{{phylo_branch}}/timetree_calibrations.log"
    benchmark: f"{LOG}/{{phylo_branch}}/benchmarks/timetree_calibrations.tsv"
    shell:
        "{PYTHON:q} {input.code:q} --tree {input.tree:q} --samples {input.samples:q} "
        "--taxonomy {input.taxonomy:q} --outdir {params.outdir:q} > {log:q} 2>&1"


rule date_busco_species_tree:
    wildcard_constraints: phylo_branch=DATING_BRANCH_PATTERN
    input:
        tree=f"{PHYLO_RUN}/species_tree.nwk",
        provenance=f"{PHYLO_RUN}/species_tree.json",
        calibrations=f"{PHYLO_RUN}/{CALIBRATION_SOURCE}/calibrations.tsv",
        calibration_reports=[f"{PHYLO_RUN}/{name}" for name in CALIBRATION_PATHS[1:]],
        code=f"{SCRIPTS}/date_phylogeny.py",
        helpers=[f"{SCRIPTS}/{name}.py" for name in
                 ["infer_phylogeny", "busco_phylogeny", "common", "sample_identity", "phylogeny_outgroup", "nwkit_runtime"]]
    output:
        tree=f"{PHYLO_RUN}/dating/species_tree.dated.nwk",
        provenance=f"{PHYLO_RUN}/dating/provenance.json",
        ages=f"{PHYLO_RUN}/dating/node_ages.tsv",
        config=f"{PHYLO_RUN}/dating/treepl.conf",
        native_log=f"{PHYLO_RUN}/dating/treepl.log"
    params: outdir=f"{PHYLO_RUN}/dating", source=CALIBRATION_SOURCE
    conda: "../envs/dating.yaml"
    threads: 1
    resources: mem_mb=4000
    log: f"{LOG}/{{phylo_branch}}/dating.log"
    benchmark: f"{LOG}/{{phylo_branch}}/benchmarks/dating.tsv"
    shell:
        "{PYTHON:q} {input.code:q} --tree {input.tree:q} --provenance {input.provenance:q} "
        "--calibrations {input.calibrations:q} --calibration-source {params.source:q} "
        "--outdir {params.outdir:q} > {log:q} 2>&1"


rule plot_dated_species_tree:
    wildcard_constraints: phylo_branch=DATING_BRANCH_PATTERN
    input:
        tree=rules.date_busco_species_tree.output.tree,
        ages=rules.date_busco_species_tree.output.ages,
        samples=dating_samples,
        traits=rules.normalize_traits.output.traits,
        taxonomy=TAXONOMY_DB,
        code=f"{SCRIPTS}/plot_dated_tree.py",
        helpers=[f"{SCRIPTS}/{name}.py" for name in
                 ["common", "infer_phylogeny", "sample_identity", "species_traits", "taxonomy_check",
                  "busco_phylogeny", "phylogeny_outgroup"]]
    output:
        pdf=f"{PHYLO_RUN}/dating/species_tree.dated.pdf",
        svg=f"{PHYLO_RUN}/dating/species_tree.dated.svg",
        png=f"{PHYLO_RUN}/dating/species_tree.dated.png",
        tips=f"{PHYLO_RUN}/dating/species_tree.dated.tips.tsv"
    params: outdir=f"{PHYLO_RUN}/dating", trait=PHY["trait"]
    conda: "../envs/analysis.yaml"
    threads: 1
    resources: mem_mb=1000
    log: f"{LOG}/{{phylo_branch}}/dating_plot.log"
    shell:
        "{PYTHON:q} {input.code:q} --tree {input.tree:q} --ages {input.ages:q} --samples {input.samples:q} "
        "--taxonomy {input.taxonomy:q} --traits {input.traits:q} "
        "--trait {params.trait:q} --outdir {params.outdir:q} > {log:q} 2>&1"
