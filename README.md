# phenoradar_prep

A Snakemake workflow that turns transcriptome assemblies and abundance estimates
into orthogroup-level TPM tables for comparative expression analysis, with an
optional KEGG Orthology (KO) annotation and expression branch.
An optional BUSCO protein phylogeny branch uses cdskit, FAMSA, trimAl, VeryFastTree, and
ASTRAL-IV/CASTLES-II, with manual or nwkit/TimeTree secondary calibrations for
optional LSD2 dating.
`phylogeny.species_sets` selects all BUSCO-filtered species, all species with an
observed phenotype, or both, retaining independent outputs for later reuse.
The optional `taxonomy_audit` target uses MonoPhy to report unexpected species-tree
placements for manual review, with group membership, run IDs and figures at each rank.
The `filter_species` target exports a subset of completed outputs using the
top-level `exclude_species` list, preserving the original analysis.
The `contrast_pairs` target selects trait-guided representatives with nwkit,
reuses the BUSCO inference rules, and exports contrast pairs and a summary tree.
`phylogeny_contrast_pairs` assigns pairs directly from full/phenotyped molecular
trees. Species exclusion recomputes those pairs inside `filtered/`, using the
completed trees without rerunning inference or the NCBI representative analysis.
Each inference run can select an outgroup automatically within its own input species set.
An optional OG alignment branch saves untrimmed protein alignments for every
observed orthogroup, retaining all gene copies for downstream site analysis.
The explicit `phenoradar_inputs` target validates completed outputs and collects
only selected PhenoRadar inputs as links, including independent KO expression
and KO-to-module/pathway maps.

```text
Metadata + BUSCO + taxonomy -> species selection -> CDS translation
  -> ODB-mapper -> orthogroup-level TPM tables, QC, and provenance
                -> optional FAMSA -> all-copy, untrimmed OG protein alignments
  -> optional KofamScan -> KO-level TPM sums, membership tables, QC, and provenance
```

Dataset inputs are grouped under `input/`; external paths or links can also be
configured. Workflow jobs read input files without modification. Datasets,
reference databases, and analysis outputs are not distributed with the repository.

## Quick start

Use Linux with Bash and Conda. Run the following commands from the repository root:

```bash
conda env create -n phenoradar-workflow -f environment.yaml
conda activate phenoradar-workflow
```

For a new dataset, create a configuration; skip the copy if you already have one:

```bash
cp config/config.yaml config/mydata.yaml
```

Edit the metadata, BUSCO, CDS, and abundance paths, optional `taxonomy.source`, and
`odb.node`. The default node `3193` is dataset-specific. See
[inputs and configuration](docs/configuration.md) for the required formats.
Missing [taxonomy snapshots](docs/references.md#taxonomy-reference), OrthoDB
references, and optional [KOfam/KEGG references](docs/kegg.md) are downloaded
and prepared automatically when needed. Existing
taxonomy snapshots are reused without updates.

Submit the full workflow as one Slurm job, replacing `YOUR_PARTITION`:

```bash
mkdir -p logs
sbatch --partition=YOUR_PARTITION \
  run_pipeline.sh --configfile config/mydata.yaml
```

All steps run locally inside that allocation, and the terminal can be closed after
submission. The default allocation is 16 CPUs and 192 GiB, sized for one ODB chunk.
Set total resources with `sbatch` options or the script's `#SBATCH` lines.
For direct execution on the current host, use the same script with a CPU and memory budget:

```bash
./run_pipeline.sh --software-deployment-method conda \
  --configfile config/mydata.yaml --cores 16 \
  --resources mem_gb=192
```

To run only species selection and manifest preparation, append `-- prepare`.
The [running guide](docs/running.md) covers smaller preparation allocations,
dry-runs, pilots, monitoring, and resuming an interrupted run.

## Results and documentation

Main tables are written to `results/<analysis>/orthogroups/expression/`; logs and
temporary work go to `logs/<analysis>/` and `work/<analysis>/`. The default analysis name is `full`;
`pilot` is a separate analysis using a selected subset of species.

| Guide | Contents |
| --- | --- |
| [Inputs and configuration](docs/configuration.md) | File formats, species selection, configuration keys, and path handling |
| [Running and resuming](docs/running.md) | Installation, Slurm and direct execution, resource budgets, pilots, and recovery |
| [Reference data](docs/references.md) | Taxonomy, OrthoDB and KEGG downloads, fixed storage, and reference updates |
| [KEGG annotation and KO expression](docs/kegg.md) | Optional KofamScan branch, frozen reference setup, ambiguity, and PhenoRadar inputs |
| [OG protein alignments](docs/alignments.md) | All OGs and copies, untrimmed FASTA with species-encoded gene IDs, and resuming per OG |
| [BUSCO phylogeny and dating](docs/phylogeny.md) | Existing BUSCO/CDS inputs, FAMSA, VeryFastTree, ASTRAL-IV, resources, and optional time trees |
| [Taxonomic review](docs/taxonomy_audit.md) | MonoPhy species-tree taxonomy review, intruder/outlier flags, run IDs and rank figures |
| [Manual species exclusion](docs/species_filter.md) | Top-level exclusion list, reusable filtered outputs, unchanged expression values/sites, and pruned trees |
| [PhenoRadar inputs](docs/phenoradar_inputs.md) | Minimal metadata, completed expression/sequence/tree inputs, optional KO groups, and validation |
| [Directory layout](docs/directory_layout.md) | Input organization, shared output paths, and migration records |
| [Contrast pairs](docs/contrast_pairs.md) | Representative analysis, pairs from full/phenotyped trees, recomputation after exclusion, membership and figures |
| [Phylogeny method review](docs/phylogeny_comparison.md) | Step-by-step comparison with GeneGalleon, scientific limitations, and improvement priorities |
| [Outputs and TPM interpretation](docs/outputs.md) | Generated files, normalization, ambiguous mappings, and QC |
| [Testing and validation](docs/development.md) | Test requirements, coverage, and validation limits |
