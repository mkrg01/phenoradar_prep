# phenoradar_prep

A Snakemake workflow that selects species by BUSCO completeness, maps transcriptome
CDS to OrthoDB orthogroups, and aggregates expression for
[PhenoRadar](https://github.com/mkrg01/phenoradar). Optional analyses provide KO
assignments, protein alignments, BUSCO species trees, and trait contrast pairs.
See [input preparation and formats](docs/inputs.md) for upstream data from
AMALGKIT and GeneGalleon.

## Quick start

Use Linux with Bash, Conda, and Singularity (or Apptainer). The launcher runs
Snakemake on the host and analysis tools in the container by default.
From the repository root:

```bash
conda env create -n phenoradar-workflow -f environment.yaml
conda activate phenoradar-workflow
cp config/config.yaml config/mydata.yaml
```

Edit `config/mydata.yaml` for your dataset using the [input guide](docs/inputs.md)
and [configuration guide](docs/configuration.md). Set `container_image` to a
release image URI or an absolute SIF path; see [container setup](docs/containers.md).
Follow [running the workflow](docs/running.md) for direct execution or Slurm
submission, adjusting resource settings for your computing environment.

## Results

Main expression tables are in `results/<run_name>/orthogroups/expression/`,
with `run_name: run001` by default. Input files are read without modification.
Logs and temporary work are stored in `logs/<run_name>/` and `work/<run_name>/`.

[Documentation](docs/index.md) covers configuration, optional analyses, output
formats, and collecting PhenoRadar inputs.
