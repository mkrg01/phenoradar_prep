# phenoradar_prep

Build cross-species RNA-seq datasets for [PhenoRadar](https://github.com/mkrg01/phenoradar).
The workflow assembles transcriptomes, quantifies gene expression, and stores
orthogroup expression in a reusable database. It can then use that database
to generate protein alignments, KEGG annotations, species trees, and trait
contrast pairs.

## Quick start

Use Linux x86-64 with Conda, Apptainer or Singularity, and Slurm or AGE. Follow the
[installation guide](docs/running.md#installation-and-normal-execution) to set
up the environment and container. Run the commands below from the repository root.

1. Prepare [sample metadata](docs/inputs.md) and edit [config/build.yaml](config/build.yaml).
   Set `reuse_from: null` for a first build.
2. Preview and submit the build. Review `reports/plan/plan.tsv` after `plan`.

   ```bash
   ./run_build.sh plan
   ./run_build.sh submit
   ```

3. After the build completes, edit [config/analysis.yaml](config/analysis.yaml)
   for your cohort, trait, and optional analyses. Use the exact build name
   printed by `submit` and choose an analysis name.

   ```bash
   ./run_analysis.sh plan --build results/BUILD
   ./run_analysis.sh submit --build results/BUILD --name ANALYSIS
   ```

The reusable database is in `results/BUILD/database/`. The
[PhenoRadar inputs](docs/phenoradar_inputs.md) are in
`results/BUILD/downstream/ANALYSIS/`. See [status and retries](docs/datasets.md#failure-and-recovery-behavior)
for saved runs.

[Documentation](docs/index.md) · [Curated datasets](datasets/README.md) ·
[Angiosperm leaf dataset](datasets/angiosperm_leaf/README.md)
