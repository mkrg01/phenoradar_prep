# phenoradar_prep

[![Snakemake](https://img.shields.io/badge/snakemake-≥9.0.0-brightgreen.svg)](https://snakemake.github.io)

Prepare [PhenoRadar](https://github.com/mkrg01/phenoradar) inputs from
[AMALGKIT](https://github.com/kfuku52/amalgkit) RNA-seq metadata.
[GeneGalleon](https://github.com/kfuku52/genegalleon) retrieves reads, assembles
transcripts, runs BUSCO, and quantifies expression. This workflow maps proteins
to OrthoDB and produces OG expression, with optional alignments, KO annotation,
and phylogenetic analyses. Completed sample data can be reused across builds.

## Requirements

Linux x86-64 with Conda, Slurm, and Apptainer/Singularity. See
[installation](docs/running.md#installation-and-normal-execution) for host setup.

## Quick start

Prepare [input/metadata.tsv](docs/inputs.md), then edit
[build.yaml](config/build.yaml) and [analysis.yaml](config/analysis.yaml).
For a first build, set `reuse_from: null`. Choose the trait and optional analyses
before submitting; the supplied analysis config enables representative trees
and contrast pairs.

Run from the repository root:

```bash
conda env create -n phenoradar_prep -f environment.yaml
conda activate phenoradar_prep
./run_pipeline.sh --prepare-container --cores 1 --resources mem_gb=4

./run_build.sh plan
./run_build.sh submit --name leaf

# After the build completes:
./run_analysis.sh submit --build results/leaf --name carnivory
```

The first submission saves inputs and settings. Resubmit the same name to retry;
use a new name when changing samples or scientific settings.

Reusable products are in `results/leaf/database/`. Analysis results, including
`phenoradar_inputs/`, are in `results/leaf/downstream/carnivory/`.

For maintained representative datasets, see the [dataset catalog](datasets/README.md)
and [angiosperm leaf update instructions](datasets/angiosperm_leaf/README.md).
Metadata acquisition, curation, review, and acceptance use `run_metadata.sh`;
database builds start separately with `run_build.sh`. Dataset metadata can stay
in `datasets/` and be referenced directly by the build config.

[Documentation](docs/index.md) · [Build and analysis guide](docs/datasets.md) ·
[Resources and retries](docs/running.md) · [Outputs](docs/outputs.md)
