# Running and resuming the workflow

[Back to README](../README.md)

Run all commands from the repository root. First configure your [inputs](configuration.md)
and provide the required [taxonomy snapshot](references.md#taxonomy-reference).

## Installation

Use Linux with Bash and Conda available. Run these commands from the repository root:

```bash
conda env create -n phenoradar-workflow -f environment.yaml
conda activate phenoradar-workflow
```

This environment provides Snakemake. For direct execution, add
`--software-deployment-method conda` to create the separate software environments
defined in `workflow/envs/` for each processing step. Inside a Slurm allocation,
the launcher enables this option automatically.

The environment definitions pin the main package versions, including Snakemake
9.8.0 and Orthologer 3.8.1. They are not complete dependency lockfiles.

For direct execution with existing software environments, set `tools.python`,
`tools.seqkit`, `tools.odb_command`, and, when applicable, `tools.odb_prefix` in
your configuration. Use `run_pipeline.sh` without
`--software-deployment-method conda`. Keep paths in your own configuration;
no machine-specific configuration is distributed with the workflow.
The installed ODB-mapper command must support OrthoDB v12.

## Execution modes

Use the same `run_pipeline.sh` script for both modes:

- `sbatch run_pipeline.sh ...` reserves one Slurm allocation for the whole workflow.
- `./run_pipeline.sh ...` runs the workflow on the current host.

The `#SBATCH` lines at the top of the script set batch defaults. Bash treats them
as comments during direct execution, so they do not reserve resources in that mode.
Snakemake uses its local executor in both modes.

Choose a target to stop at a particular stage. Put targets after `--`, with all
options before it, for example `--configfile config/mydata.yaml -- mapping`.
This also keeps `--configfile` from interpreting the target as another filename.

| Target | Result |
| --- | --- |
| `all` (default) | Complete workflow through the combined TPM and QC tables |
| `prepare` | Selected metadata, provenance, and sample/chunk manifests |
| `references` | Prepared OrthoDB reference snapshot |
| `proteins` | Translated CDS for the selected species |
| `mapping` | ODB chunk results and the merged gene-to-orthogroup index |
| `kegg` | KO annotations and original-TPM sums; independent of ODB mapping, after [KEGG reference setup](kegg.md) |

Prerequisite steps are included automatically. See [reference setup](references.md)
and [output formats](outputs.md) for details.

## Slurm: run the workflow in one allocation

`sbatch run_pipeline.sh` reserves one Slurm allocation and runs Snakemake with its
**local executor** inside that allocation. All processing steps run as ordinary
processes on the allocated node. They share its CPU and memory budget and do not
submit additional Slurm jobs. You can close the terminal after `sbatch` returns.

From the repository root, activate the environment and create the log directory
**before** submission. Replace `YOUR_PARTITION` with your cluster's partition:

```bash
conda activate phenoradar-workflow
mkdir -p logs

# Run the full workflow: 32 CPUs, 256 GiB, up to 21 days by default.
sbatch --partition=YOUR_PARTITION \
  run_pipeline.sh --configfile config/mydata.yaml
```

To run only preparation, request a smaller allocation:

```bash
sbatch --partition=YOUR_PARTITION --cpus-per-task=2 --mem=8G --time=01:00:00 \
  run_pipeline.sh --configfile config/mydata.yaml -- prepare
```

Preparation is also included automatically in the full workflow. Do not start
both commands concurrently for the same outputs. Arguments before the script
name configure the Slurm allocation; arguments after it specify workflow options
and targets. The activated environment is exported to the job, and Conda manages
the processing environments automatically.

| Setting | Meaning in this mode |
| --- | --- |
| `#SBATCH --cpus-per-task=32` | CPU budget shared by all concurrently running steps; passed to Snakemake as `--cores` |
| `#SBATCH --mem=256G` | Total allocation memory, including Snakemake and its processes |
| `#SBATCH --time=21-00:00:00` | Time limit for the entire workflow after the allocation starts |
| Rule `threads`, including `odb.threads` | Threads used by a step, capped by the allocation's CPU budget |
| `odb.mem_gb` | Memory estimate for one mapping step, in GB; converted to Snakemake's `mem_mb` scheduling resource |
| Internal `odb_slots=1` | At most one ODB reference/mapping step at a time; set automatically by the launcher |

The script derives the memory budget from Slurm's allocation and leaves 4 GB
for Snakemake and environment management. It displays the budget in GB.
Declared rule memory values are scheduling estimates; Slurm enforces the overall
allocation limit. Use one node, one task, and an explicit finite memory request.
Both `--mem` and `--mem-per-cpu` are supported. Change the total budget with
`sbatch` options, for example:

```bash
sbatch --partition=YOUR_PARTITION --cpus-per-task=32 --mem=192G --time=21-00:00:00 \
  run_pipeline.sh --configfile config/mydata.yaml
```

These defaults are initial allowances, not measured requirements. The allocation
must fit the largest processing step and last long enough for the whole workflow,
including all sequential chunks. Reserving resources once avoids queue waits
between steps, but keeps those resources reserved during lighter steps as well.
A larger or longer allocation may also wait longer before starting.

The entry point enforces local execution and the allocation budgets. Specify
the partition and any required account with `sbatch --partition=...` and
`sbatch --account=...`. The compute node needs access to the configured files and
any required package/download services. Choose paths consistent with the
workflow's storage checks.

`sbatch` prints a single job ID. Monitor or stop the whole workflow with:

```bash
squeue -u "$USER"
tail -f logs/pipeline-JOB_ID.log
scancel JOB_ID
```

Replace `JOB_ID` with the printed ID. Stopping this allocation stops its processing
steps. Resubmit the same command to resume missing or incomplete outputs. The
script requests an interrupt shortly before the time limit for orderly cleanup;
it preserves Snakemake's exit status as the batch job's exit status.

## Direct execution

Outside Slurm, specify an appropriate CPU and memory budget:

```bash
./run_pipeline.sh --software-deployment-method conda \
  --configfile config/mydata.yaml --cores 24 \
  --resources mem_gb=128
```

This runs processing on the current host. The memory value is a scheduling budget,
so leave additional physical memory for Snakemake and other processes. For shared
clusters, submit the same script with `sbatch` to reserve resources and keep running
after terminal disconnection.

The launcher sets `odb_slots=1` automatically, including when you supply
`--resources mem_gb=...`. You do not need to repeat it in commands. This limits
ODB concurrency within Snakemake; it does not submit or reserve Slurm jobs.

When `SLURM_JOB_ID` is present, including direct execution inside an existing
Slurm allocation, the script uses that allocation's CPU and memory budgets and
enables Conda automatically. Run from the repository root in that case; set the
total resources when requesting the allocation.

## Memory units

Use positive whole numbers for `odb.mem_gb` in the configuration and
`--resources mem_gb=...` when launching directly. Both use decimal GB:
`128` means 128 GB, equivalent to 128000 MB. `odb.mem_gb` is the estimate for one
mapping step; the command-line resource is the total budget for concurrent steps.

The launcher and workflow convert these values to Snakemake's standard `mem_mb`
resource internally. Snakemake's own job listings can therefore still show MB.
Slurm's `--mem=256G` uses GiB (1 GiB = 1024 cubed bytes). The launcher accounts for
that difference when calculating the budget from the allocation; no manual
conversion is needed when submitting with `sbatch`.

## Prepare and inspect

To select species and create the sample and chunk manifests without downloading
OrthoDB or running ODB-mapper, use the `prepare` batch command above. Species
selection is a checkpoint: a dry-run before preparation may not display all
downstream jobs. After `prepare` finishes, inspect the complete plan with:

```bash
./run_pipeline.sh --software-deployment-method conda \
  --configfile config/mydata.yaml --cores 24 --dry-run
```

## Pilot run

After preparation, create `config/pilot_species.txt` with a small set of eligible
species from `results/<analysis>/metadata/species_high_busco.txt`. For example,
with the default analysis name `full`, select the first three species and review
or edit the list:

```bash
head -n 3 results/full/metadata/species_high_busco.txt > config/pilot_species.txt
```

The generic `config/pilot.yaml` reads that list, uses two species per ODB chunk,
retains work directories, and writes to `results/pilot/`. A missing list causes
an input error. Run the pilot in one allocation with:

```bash
sbatch --partition=YOUR_PARTITION --cpus-per-task=8 --mem=80G \
  run_pipeline.sh --configfile config/mydata.yaml config/pilot.yaml
```

## Re-running and recovery

Run the same command again to process missing or outdated outputs. Snakemake
tracks declared inputs, parameters, code, and software environment definitions.
Protein FASTA files are individually declared as inputs to ODB jobs. Changing only
an abundance TSV recalculates that run's TPM and the combined tables.

```bash
# Resume after an interruption.
sbatch --partition=YOUR_PARTITION \
  run_pipeline.sh --configfile config/mydata.yaml

# Explicitly rerun mapping after manually changing installed ODB software.
sbatch --partition=YOUR_PARTITION \
  run_pipeline.sh --configfile config/mydata.yaml --forcerun odb_map
```

Incomplete ODB work is retained under
`work/<analysis>/odb/<chunk>/<fingerprint>/`. The fingerprint includes input FASTA
contents, the reference record, mapping settings, software records, and worker
code. Matching incomplete work can be resumed; changed inputs or settings use a
different work directory. A completed result is not skipped merely because a
success marker exists.

After a successful job, native results and logs are copied out before that job's
work directory is removed. Set `odb.keep_work: true` to retain it, as the pilot
configuration does. Work from other fingerprints and older outputs are not
automatically deleted. Final tables use only the species and runs in the current
sample manifest, even when older per-run outputs remain on disk.
