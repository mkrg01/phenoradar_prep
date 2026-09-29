# phenoradar_prep

[![Snakemake](https://img.shields.io/badge/snakemake-≥9.0.0-brightgreen.svg)](https://snakemake.github.io)

Prepare [PhenoRadar](https://github.com/mkrg01/phenoradar) inputs from manually
curated [AMALGKIT](https://github.com/kfuku52/amalgkit) RNA-seq metadata.
[GeneGalleon](https://github.com/kfuku52/genegalleon) handles assembly, BUSCO, and
quantification; Snakemake maps CDS to OrthoDB and runs downstream analyses.
Adding or removing species reuses completed species products and mappings.

## Requirements

- [Workflow environment](environment.yaml), including Snakemake and its Slurm executor
- Apptainer / Singularity
- Slurm

## Quick start

Prepare `input/metadata.tsv` with one run per independent sample, from NCBI or local FASTQ
files. Edit [build.yaml](config/build.yaml) and [analysis.yaml](config/analysis.yaml)
directly. See the [initial setup](docs/running.md#installation-and-normal-execution) and
[build and analysis guide](docs/datasets.md).
Use the `name` from `config/build.yaml` in the build paths below.

```bash
# Cache the workflow container: first run, new release, or cleared image cache.
./run_pipeline.sh --prepare-container --cores 1 --resources mem_gb=4

# Preview reusable products, missing work, and conflicts.
./run_build.sh plan

# Save build inputs and settings in builds/<name>.
./run_build.sh prepare

# Submit assembly, BUSCO QC, expression quantification, and OrthoDB mapping to Slurm.
./run_build.sh submit --build builds/angiosperm_leaf_20260925 --until mapping

# After the build finishes, select samples and save analysis inputs and settings.
./run_analysis.sh prepare --build builds/angiosperm_leaf_20260925 --name analysis001

# Submit expression aggregation and enabled analyses; collect PhenoRadar inputs.
./run_analysis.sh submit --analysis analyses/analysis001
```

Builds can stop after assembly, BUSCO, or quantification for inspection.
Its `products/` directory can be copied to another project for reuse.
Use a new build name after editing metadata; removed species leave the new outputs.

Final inputs are written to `results/<analysis>/phenoradar_inputs/`.
See the [documentation](docs/index.md) for configuration and output formats.
