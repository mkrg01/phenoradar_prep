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

Edit [config/config.yaml](config/config.yaml) directly using the
[input guide](docs/inputs.md) and [configuration guide](docs/configuration.md).
For Slurm, adjust the `#SBATCH` settings in `run_pipeline.sh` for your cluster,
then submit from the repository root:

```bash
sbatch run_pipeline.sh
```

The workflow reads `config/config.yaml` automatically and uses the CPUs and
memory allocated by Slurm. See [running the workflow](docs/running.md) for direct
execution, individual targets, and resource settings.

## Results

Main expression tables are in `results/<run_name>/orthogroups/expression/`,
with `run_name: run001` by default.
Slurm writes standard output to `pipeline-<job_id>.out` and standard error to
`pipeline-<job_id>.err` in the repository root. Per-step logs and temporary work
are stored in `logs/<run_name>/` and `work/<run_name>/`.

[Documentation](docs/index.md) covers configuration, optional analyses, output
formats, and collecting PhenoRadar inputs.
