# phenoradar_prep

[![Snakemake](https://img.shields.io/badge/snakemake-≥9.0.0-brightgreen.svg)](https://snakemake.github.io)

A Snakemake workflow that selects species by BUSCO completeness, maps transcriptome
CDS to OrthoDB orthogroups, and aggregates expression for
[PhenoRadar](https://github.com/mkrg01/phenoradar). Optional analyses provide KO
assignments, protein alignments, BUSCO species trees, and trait contrast pairs.
See [input preparation and formats](docs/inputs.md) for upstream data from
[AMALGKIT](https://github.com/kfuku52/amalgkit) and [GeneGalleon](https://github.com/kfuku52/genegalleon).

## Requirements

Before running the workflow, make sure the following software is installed:

- [Snakemake ≥ 9.0.0](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
- [Apptainer (Singularity)](https://apptainer.org/docs/admin/main/installation.html)

## Quick start

From the repository root:

```bash
cp config/config.yaml config/mydata.yaml
```

Edit `config/mydata.yaml` for your dataset using the [input guide](docs/inputs.md)
and [configuration guide](docs/configuration.md).
Follow [running the workflow](docs/running.md) for direct execution or Slurm
submission, adjusting resource settings for your computing environment.

## Results

Main expression tables are in `results/<run_name>/orthogroups/expression/`,
with `run_name: run001` by default.
Logs and temporary work are stored in `logs/<run_name>/` and `work/<run_name>/`.

[Documentation](docs/index.md) covers configuration, optional analyses, output
formats, and collecting PhenoRadar inputs.
