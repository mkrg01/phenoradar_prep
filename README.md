# phenoradar_prep

[![Snakemake](https://img.shields.io/badge/snakemake-≥9.0.0-brightgreen.svg)](https://snakemake.github.io)

Prepare [PhenoRadar](https://github.com/mkrg01/phenoradar) inputs from [AMALGKIT](https://github.com/kfuku52/amalgkit) RNA-seq metadata.
[GeneGalleon](https://github.com/kfuku52/genegalleon) handles assembly, BUSCO, and
quantification; Snakemake maps CDS to OrthoDB and runs downstream analyses.
Adding or removing samples reuses completed sample data, mappings, and OG expression.

## Requirements

- Conda
- Apptainer / Singularity
- Slurm

## Quick start

Prepare `input/metadata.tsv` using [AMALGKIT](https://github.com/kfuku52/amalgkit).
Edit [build.yaml](config/build.yaml) and [analysis.yaml](config/analysis.yaml)
directly. See the [initial setup](docs/running.md#installation-and-normal-execution) and
[build and analysis guide](docs/datasets.md).
Use the `name` from `config/build.yaml` in the build paths below.

```bash
# From the repository root: create once, activate each session.
conda env create -n phenoradar_prep -f environment.yaml
conda activate phenoradar_prep

# Cache the workflow container: first run, new release, or cleared image cache.
./run_pipeline.sh --prepare-container --cores 1 --resources mem_gb=4

# Save the build plan to results/build_plan.tsv.
mkdir -p results
./run_build.sh plan > results/build_plan.tsv

# Save inputs and settings, then submit the build through OG expression.
./run_build.sh submit

# After the build finishes, save conditions and submit downstream analyses.
./run_analysis.sh submit --build results/angiosperm_leaf_20260925 --name analysis001
```

Under `results/<name>/`, `database/` holds reusable sample data,
`downstream/<condition>/` holds results for selected samples (including
`phenoradar_inputs/`), and `work/` and `logs/` retain computation files.

Builds can stop after assembly, BUSCO, or quantification. After mapping, OG
expression, and validation finish, `database/` can be copied to another project.
Submitting the same name reuses saved conditions; use a new name after changing
inputs or settings. See the [documentation](docs/index.md) for configuration and output formats.
