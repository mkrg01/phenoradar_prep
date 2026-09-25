# phenoradar_prep

[![Snakemake](https://img.shields.io/badge/snakemake-≥9.0.0-brightgreen.svg)](https://snakemake.github.io)

Prepare [PhenoRadar](https://github.com/mkrg01/phenoradar) inputs from manually
curated [AMALGKIT](https://github.com/kfuku52/amalgkit) RNA-seq metadata.
[GeneGalleon](https://github.com/kfuku52/genegalleon) handles assembly, BUSCO, and
quantification; Snakemake maps CDS to OrthoDB and runs downstream analyses.
Adding or removing species reuses completed species products and mappings.

## Requirements

- [Workflow environment](environment.yaml), including Snakemake and its Slurm executor
- [Apptainer / Singularity](docs/running.md#installation-and-normal-execution)
- Slurm

## Quick start

Prepare `input/metadata.tsv` with one run per species, from NCBI or local FASTQ
files. Edit [build.yaml](config/build.yaml) and [analysis.yaml](config/analysis.yaml)
directly. Set `name` in build settings (e.g. `angiosperm_leaf_20260925`);
`prepare --name` overrides it. See the [build and analysis guide](docs/datasets.md).

```bash
./run_build.sh plan
./run_build.sh prepare
./run_build.sh submit --build builds/angiosperm_leaf_20260925 --until mapping

# After the build finishes:
./run_analysis.sh prepare --build builds/angiosperm_leaf_20260925 --name analysis001
./run_analysis.sh submit --analysis analyses/analysis001
```

Builds can stop after assembly, BUSCO, or quantification for inspection. A completed
build supports multiple analyses with different species selections and traits.
Its `products/` directory can be copied to another project for reuse.
Use a new build name after editing metadata; removed species leave the new outputs.

Final inputs are written to `results/<analysis>/phenoradar_inputs/`.
See the [documentation](docs/index.md) for configuration and output formats.
