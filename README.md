# phenoradar_prep

[![Snakemake](https://img.shields.io/badge/snakemake-≥9.0.0-brightgreen.svg)](https://snakemake.github.io)

A Snakemake workflow that selects species by BUSCO completeness, maps transcriptome
CDS to OrthoDB orthogroups, and aggregates expression for
[PhenoRadar](https://github.com/mkrg01/phenoradar). Optional analyses provide KO
assignments, protein alignments, BUSCO species trees, and trait contrast pairs.
See [input preparation and formats](docs/inputs.md) for upstream data from
[AMALGKIT](https://github.com/kfuku52/amalgkit) and [GeneGalleon](https://github.com/kfuku52/genegalleon).

To add RNA-seq species incrementally, use the [dataset interface](docs/datasets.md):
manual metadata, reuse of completed products, staged GeneGalleon Slurm arrays,
and ODB mapping of missing species only. Removing metadata rows removes species
from the next dataset output while retaining reusable products. Missing pinned
GeneGalleon source and SIF dependencies are acquired automatically before new work.

## Requirements

Before running the workflow, make sure the following software is installed:

- [Snakemake ≥ 9.0.0](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
- [Apptainer (Singularity)](https://apptainer.org/docs/admin/main/installation.html)

## Quick start

Place prepared files in `input/` following the [input guide](docs/inputs.md),
then edit [config/config.yaml](config/config.yaml) using the
[configuration guide](docs/configuration.md).
For Slurm, adjust the `#SBATCH` settings in `run_pipeline.sh` for your cluster,
then submit from the repository root:

```bash
sbatch run_pipeline.sh
```

Choose species trees with `phylogeny.trees`: `[]` disables inference, `[all]`
uses every selected species, `[phenotyped]` uses known-trait species, and
`[representatives]` compresses known-trait species before inference. Pair
selection is controlled by `phylogeny.contrast_pairs.enabled` and uses those
same trees. See [configuration](docs/configuration.md#choosing-species-trees).

The workflow reads `config/config.yaml` automatically and uses the CPUs and
memory allocated by Slurm. See [running the workflow](docs/running.md) for direct
execution, individual targets, and resource settings.

## Results

Main expression tables are in `results/<run_name>/orthogroups/expression/`,
with `run_name: run001` by default.
Per-step logs and temporary work
are stored in `logs/<run_name>/` and `work/<run_name>/`.

[Documentation](docs/index.md) covers configuration, optional analyses, output
formats, and collecting PhenoRadar inputs.
