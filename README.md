# phenoradar_prep

A Snakemake workflow that turns transcriptome assemblies and abundance estimates
into orthogroup-level TPM tables for comparative expression analysis.

```text
Metadata + BUSCO + taxonomy -> species selection -> CDS translation
  -> ODB-mapper -> orthogroup-level TPM tables, QC, and provenance
```

Input files are read without modification. Datasets, reference databases, and
analysis outputs are not distributed with the repository.

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

Edit the metadata, BUSCO, CDS, and abundance paths, `taxonomy.database`, and
`odb.node`. The default node `3193` is dataset-specific. See
[inputs and configuration](docs/configuration.md) for the required formats.
Provide a [taxonomy snapshot](docs/references.md#taxonomy-reference) before running;
OrthoDB is downloaded and prepared automatically when needed.

Submit the full workflow as one Slurm job, replacing `YOUR_PARTITION`:

```bash
mkdir -p logs
sbatch --partition=YOUR_PARTITION \
  run_pipeline.sh --configfile config/mydata.yaml
```

All steps run locally inside that allocation, and the terminal can be closed after
submission. Set total resources with `sbatch` options or the script's `#SBATCH` lines.
For direct execution on the current host, use the same script with a CPU and memory budget:

```bash
./run_pipeline.sh --software-deployment-method conda \
  --configfile config/mydata.yaml --cores 24 \
  --resources mem_gb=128
```

To run only species selection and manifest preparation, append `-- prepare`.
The [running guide](docs/running.md) covers smaller preparation allocations,
dry-runs, pilots, monitoring, and resuming an interrupted run.

## Results and documentation

Main tables are written to `results/<analysis>/tpm/`; logs and temporary work go to
`logs/<analysis>/` and `work/<analysis>/`. The default analysis name is `full`;
`pilot` is a separate analysis using a selected subset of species.

| Guide | Contents |
| --- | --- |
| [Inputs and configuration](docs/configuration.md) | File formats, species selection, configuration keys, and path handling |
| [Running and resuming](docs/running.md) | Installation, Slurm and direct execution, resource budgets, pilots, and recovery |
| [Reference data](docs/references.md) | Taxonomy setup, OrthoDB downloads, storage, and reference updates |
| [Outputs and TPM interpretation](docs/outputs.md) | Generated files, normalization, ambiguous mappings, and QC |
| [Testing and validation](docs/development.md) | Test requirements, coverage, and validation limits |
