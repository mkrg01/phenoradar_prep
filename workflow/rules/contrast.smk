rule prepare_ncbi_guide:
    wildcard_constraints: guide_branch="rooting|phylogeny_phenotyped/rooting"
    input:
        samples=lambda wc: (checkpoints.select_phenotyped_species.get().output.samples
                            if wc.guide_branch == "phylogeny_phenotyped/rooting"
                            else checkpoints.select_metadata.get().output.samples),
        taxonomy=TAXONOMY_DB,
        code=f"{SCRIPTS}/phylogeny_root.py",
        helpers=f"{SCRIPTS}/common.py"
    output:
        tree=f"{OUT}/{{guide_branch}}/ncbi_tree.nwk",
        taxids=f"{OUT}/{{guide_branch}}/taxids.tsv"
    conda: "../envs/timetree.yaml"
    resources: mem_mb=4000
    log: f"{LOG}/{{guide_branch}}/ncbi.log"
    shell:
        "{PYTHON:q} {input.code:q} ncbi_tree --samples {input.samples:q} --taxonomy-db {input.taxonomy:q} "
        "--output {output.tree:q} --taxids {output.taxids:q} > {log:q} 2>&1"


rule prepare_phylogeny_outgroup:
    input:
        samples=phylogeny_samples,
        taxonomy=TAXONOMY_DB,
        metadata=f"{META}/metadata_high_busco.tsv",
        tree=phylogeny_root_guide,
        code=f"{SCRIPTS}/phylogeny_root.py",
        helpers=f"{SCRIPTS}/common.py"
    output:
        outgroup=f"{PHYLO_RUN}/rooting/outgroup.txt",
        qc=f"{PHYLO_RUN}/rooting/outgroup.json"
    params:
        outgroup=PHY["outgroup"],
        tree_flag="--tree" if PHY["outgroup"] == "auto" else ""
    conda: "../envs/timetree.yaml"
    resources: mem_mb=4000
    log: f"{LOG}/{{phylo_branch}}/rooting/outgroup.log"
    shell:
        "{PYTHON:q} {input.code:q} prepare_root --samples {input.samples:q} --metadata {input.metadata:q} "
        "--output {output.outgroup:q} --qc {output.qc:q} --outgroup {params.outgroup:q} --taxonomy-db {input.taxonomy:q} "
        "{params.tree_flag} {input.tree:q} > {log:q} 2>&1"


checkpoint select_contrast_representatives:
    input:
        samples=f"{META}/samples.tsv",
        metadata=f"{META}/metadata_high_busco.tsv",
        traits=config["inputs"]["species_trait"],
        tree=f"{ROOTING}/ncbi_tree.nwk",
        code=f"{SCRIPTS}/contrast_pairs.py",
        helpers=[f"{SCRIPTS}/common.py", f"{SCRIPTS}/species_traits.py", f"{SCRIPTS}/phylogeny_root.py"]
    output:
        samples=f"{CONTRAST}/selection/samples.tsv",
        traits=f"{CONTRAST}/selection/traits.tsv",
        tree=f"{CONTRAST}/selection/ncbi_skim.nwk",
        all=f"{CONTRAST}/selection/ncbi_skim.all.tsv",
        sampled=f"{CONTRAST}/selection/ncbi_skim.sampled.tsv",
        qc=f"{CONTRAST}/selection/selection.json"
    params:
        outdir=f"{CONTRAST}/selection", trait=config["contrast"]["trait"], seed=PHY["seed"]
    conda: "../envs/timetree.yaml"
    resources: mem_mb=4000
    log: f"{LOG}/contrast/selection.log"
    shell:
        "{PYTHON:q} {input.code:q} prepare --samples {input.samples:q} --metadata {input.metadata:q} "
        "--traits {input.traits:q} --tree {input.tree:q} "
        "--trait {params.trait:q} --seed {params.seed} --outdir {params.outdir:q} > {log:q} 2>&1"


rule identify_contrast_pairs:
    input:
        tree=f"{CONTRAST}/phylogeny/species_tree.nwk",
        outgroup=f"{CONTRAST}/phylogeny/rooting/outgroup.txt",
        selection=rules.select_contrast_representatives.output,
        code=f"{SCRIPTS}/contrast_pairs.py",
        helpers=[f"{SCRIPTS}/common.py", f"{SCRIPTS}/species_traits.py", f"{SCRIPTS}/phylogeny_root.py"]
    output:
        tree=f"{CONTRAST}/summary_tree.nwk",
        all=f"{CONTRAST}/summary_tree.all.tsv",
        sampled=f"{CONTRAST}/summary_tree.sampled.tsv",
        contrastive=f"{CONTRAST}/contrastive.nwk",
        contrast_all=f"{CONTRAST}/contrastive.all.tsv",
        contrast_sampled=f"{CONTRAST}/contrastive.sampled.tsv",
        pairs=f"{CONTRAST}/contrast_pairs.tsv",
        metadata=f"{CONTRAST}/species_metadata.tsv",
        qc=f"{CONTRAST}/summary.json"
    params:
        selection=f"{CONTRAST}/selection", outdir=CONTRAST, seed=PHY["seed"]
    conda: "../envs/timetree.yaml"
    resources: mem_mb=4000
    log: f"{LOG}/contrast/pairs.log"
    shell:
        "{PYTHON:q} {input.code:q} summarize --tree {input.tree:q} --selection-dir {params.selection:q} "
        "--outgroup-file {input.outgroup:q} --outdir {params.outdir:q} --seed {params.seed} > {log:q} 2>&1"


rule plot_contrast_tree:
    wildcard_constraints: contrast_branch="contrast|phylogeny/contrast|phylogeny_phenotyped/contrast"
    input:
        tree=f"{OUT}/{{contrast_branch}}/summary_tree.nwk",
        metadata=f"{OUT}/{{contrast_branch}}/species_metadata.tsv",
        summary=f"{OUT}/{{contrast_branch}}/summary.json",
        code=f"{SCRIPTS}/plot_contrast_tree.py",
        helpers=f"{SCRIPTS}/common.py"
    output:
        pdf=f"{OUT}/{{contrast_branch}}/summary_tree.pdf", svg=f"{OUT}/{{contrast_branch}}/summary_tree.svg"
    params: outdir=f"{OUT}/{{contrast_branch}}"
    conda: "../envs/timetree.yaml"
    resources: mem_mb=4000
    log: f"{LOG}/{{contrast_branch}}/plot.log"
    shell:
        "{PYTHON:q} {input.code:q} --tree {input.tree:q} --metadata {input.metadata:q} "
        "--summary {input.summary:q} --outdir {params.outdir:q} > {log:q} 2>&1"


rule identify_phylogeny_contrast_pairs:
    wildcard_constraints: phylo_branch="phylogeny|phylogeny_phenotyped"
    input:
        tree=f"{PHYLO_RUN}/species_tree.nwk",
        tree_qc=f"{PHYLO_RUN}/species_tree.json",
        samples=lambda wc: f"{PHENOTYPED}/selection/samples.tsv" if wc.phylo_branch == "phylogeny_phenotyped" else f"{META}/samples.tsv",
        metadata=f"{META}/metadata_high_busco.tsv",
        traits=config["inputs"]["species_trait"],
        code=f"{SCRIPTS}/contrast_pairs.py",
        helpers=[f"{SCRIPTS}/common.py", f"{SCRIPTS}/species_traits.py", f"{SCRIPTS}/phylogeny_root.py"]
    output:
        observed=f"{PHYLO_RUN}/contrast/observed_tree.nwk",
        tree=f"{PHYLO_RUN}/contrast/summary_tree.nwk",
        all=f"{PHYLO_RUN}/contrast/summary_tree.all.tsv",
        sampled=f"{PHYLO_RUN}/contrast/summary_tree.sampled.tsv",
        contrastive=f"{PHYLO_RUN}/contrast/contrastive.nwk",
        contrast_all=f"{PHYLO_RUN}/contrast/contrastive.all.tsv",
        contrast_sampled=f"{PHYLO_RUN}/contrast/contrastive.sampled.tsv",
        pairs=f"{PHYLO_RUN}/contrast/contrast_pairs.tsv",
        metadata=f"{PHYLO_RUN}/contrast/species_metadata.tsv",
        qc=f"{PHYLO_RUN}/contrast/summary.json"
    params: outdir=f"{PHYLO_RUN}/contrast", trait=config["contrast"]["trait"], seed=PHY["seed"]
    conda: "../envs/timetree.yaml"
    resources: mem_mb=4000
    log: f"{LOG}/{{phylo_branch}}/contrast/pairs.log"
    shell:
        "{PYTHON:q} {input.code:q} from_tree --tree {input.tree:q} --tree-qc {input.tree_qc:q} "
        "--samples {input.samples:q} --metadata {input.metadata:q} --traits {input.traits:q} "
        "--trait {params.trait:q} --seed {params.seed} --outdir {params.outdir:q} > {log:q} 2>&1"
