# phenoradar_prep

[![Snakemake](https://img.shields.io/badge/snakemake-≥9.0.0-brightgreen.svg)](https://snakemake.github.io)

Prepare [PhenoRadar](https://github.com/mkrg01/phenoradar) inputs from
[AMALGKIT](https://github.com/kfuku52/amalgkit) RNA-seq metadata using
[GeneGalleon](https://github.com/kfuku52/genegalleon). Build reusable sample
databases, then run OG expression and optional alignments, KO annotation,
and phylogenetic analyses.

## Requirements

Linux x86-64 with Conda, Slurm, and Apptainer/Singularity. See
[installation](docs/running.md#installation-and-normal-execution) for host setup.

## Quick start

Prepare [input/metadata.tsv](docs/inputs.md), then edit
[build.yaml](config/build.yaml) and [analysis.yaml](config/analysis.yaml).
For a first build, set `reuse_from: null`; review the trait and optional analyses.

Run from the repository root. Replace `{build_name}` with the build directory
name reported by `submit`, and choose `{analysis_name}`:

```bash
conda env create -n phenoradar_prep -f environment.yaml
conda activate phenoradar_prep
./run_pipeline.sh --prepare-container --cores 1 --resources mem_gb=4
./run_build.sh fetch-software

./run_build.sh plan
./run_build.sh submit

# After the build completes:
./run_analysis.sh submit --build results/{build_name} --name {analysis_name}
```

Build names include UTC timestamps by default. For retries, use the
[saved run path](docs/datasets.md#failure-and-recovery-behavior) via `--build`
or `--analysis`.

Reusable products are in `results/{build_name}/database/`. Analysis results,
including `phenoradar_inputs/`, are in
`results/{build_name}/downstream/{analysis_name}/`.

For curated datasets and metadata updates, see the [dataset catalog](datasets/README.md)
and [angiosperm leaf guide](datasets/angiosperm_leaf/README.md).

[Documentation](docs/index.md) · [Build and analysis guide](docs/datasets.md) ·
[Resources and retries](docs/running.md) · [Outputs](docs/outputs.md)
