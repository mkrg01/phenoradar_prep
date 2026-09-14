# phenoradar_prep

A Snakemake workflow for preparing comparative expression data from transcriptome
CDS assemblies and abundance estimates. It selects species by BUSCO completeness,
maps genes to OrthoDB orthogroups, and produces orthogroup TPM tables. Optional
analyses provide KO expression, protein alignments, BUSCO species trees, and
trait contrast pairs. Completed results can be collected for PhenoRadar.

## Quick start

Use Linux, Bash, and Conda. From the repository root:

```bash
conda env create -n phenoradar-workflow -f environment.yaml
conda activate phenoradar-workflow
cp config/config.yaml config/mydata.yaml
```

Edit `config/mydata.yaml` for your [input files](docs/inputs.md) and OrthoDB node.
The default `odb.node: 3193` is specific to the example dataset. Input data and
reference databases are not included; missing references are prepared on first use.

Submit to Slurm, replacing `YOUR_PARTITION` with your cluster's partition:

```bash
mkdir -p logs
sbatch --partition=YOUR_PARTITION \
  run_pipeline.sh --configfile config/mydata.yaml
```

The workflow runs within one allocation; the default request is 16 CPUs and
192 GiB. See [running the workflow](docs/running.md) for preparation, pilots,
direct execution, resource budgets, and resuming a run.

## Results

Main expression tables are in `results/<analysis>/orthogroups/expression/`,
with `analysis: full` by default. Input files are read without modification.
Logs and temporary work are stored in `logs/<analysis>/` and `work/<analysis>/`.

[Documentation](docs/index.md) covers configuration, optional analyses, output
formats, and collecting PhenoRadar inputs.
