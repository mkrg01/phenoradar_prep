# Running the workflow

[Documentation](index.md) · [Configuration](configuration.md)

## Installation and normal execution

Use Slurm on Linux x86-64/Bash with the [requirements](../README.md#requirements).
Compute nodes need host Snakemake and `singularity` on `PATH`; Apptainer must
provide its `singularity` compatibility command.

From the repository root, create the host environment once, then activate it in
each session before running or submitting workflow commands:

```bash
conda env create -n phenoradar_prep -f environment.yaml
conda activate phenoradar_prep
```

The launchers use Python and Snakemake from the active environment; they do not
activate it automatically. `run_pipeline.sh` uses `snakemake` on `PATH` unless
`SNAKEMAKE_BIN` is set. Analysis tools run inside the workflow container.
Apptainer/Singularity and Slurm must also be available on the host.

Use a published release and prepare the workflow image matching [VERSION](../VERSION)
before the first mapping or analysis job:

```bash
./run_pipeline.sh --prepare-container --cores 1 --resources mem_gb=4
```

This step is required with Snakemake 9.8: it can check Conda inside the image
before its normal automatic pull. Repeat after changing releases or clearing the
image cache; otherwise the cached image is reused. GeneGalleon's separate
[software and image](datasets.md#pinned-genegalleon-source-and-sif) are fetched automatically.

Follow the [build/analysis commands](datasets.md) for normal execution.
Keep imported products inside the repository, which is mounted automatically.

`plan` previews locally. The first `submit` saves inputs and settings, validates
them, submits Slurm jobs, and returns without waiting for completion. Submitting
the same name uses saved conditions. `submit --dry-run` saves those conditions
and previews job scripts without submitting; it does not check the full Snakemake
DAG. The optional `prepare` command saves inputs and settings without job scripts.
Use `status`, `squeue`, and the prepared run's `jobs/logs/` to inspect progress.

Assembly, BUSCO, and quantification use sample arrays with `afterok` dependencies.
Mapping and analysis use a controller with separate Snakemake worker jobs.
Each array task, controller, and worker has its own resource and time limits.

## Targets

Use `run_analysis.sh submit --analysis results/<build>/downstream/<analysis> --target <target>`.
Targets include missing prerequisites within the analysis; they never rebuild
upstream species products. Edit branch settings before the first submission, including `--dry-run`.

| Analysis target | Work requested |
| --- | --- |
| `all` (default) | Merge saved OG expression/QC and run enabled branches, then collect PhenoRadar inputs |
| `alignments` | [All-copy OG alignments](alignments.md) |
| `kegg` | [KO annotation and original-TPM sums](kegg.md) |
| `phylogeny_prepare` | BUSCO input audit, outgroup, and marker plan |
| `phylogeny` | Infer trees selected by `phylogeny.trees` |
| `phylogeny_calibrations` | Trees and [TimeTree calibrations](dating.md#timetree-calibrations), without dating |
| `timetree` | [LSD2 dating](dating.md) |
| `taxonomy_check` | [MonoPhy review](taxonomy_check.md) |
| `contrast_pairs` | [Trait pairs](contrast_pairs.md) on selected trees |
| `phenoradar_inputs` | [Collect completed results](phenoradar_inputs.md), without starting producers |

Phylogeny postprocessing targets require their enabled flag. Explicit `alignments`
and `kegg` targets work without enabling their inclusion in `all`.
For legacy bundles without saved OG expression, analysis computes it from the
build's quantification and mappings; see [legacy outputs](outputs.md#directory-layout).
Build endpoints are documented [separately](datasets.md#build-a-reusable-database).

## Resource budgets

Edit the relevant config's `slurm` section:

| Setting | Controls |
| --- | --- |
| `partition` | Slurm partition (`--partition`); null uses the cluster default |
| `array_size` (build only) | Sample index slots per array batch |
| `total_limits.jobs` | Running sample array tasks, or queued/running Snakemake workers; excludes the controller; null adds no job-count cap (default) |
| `total_limits.cpus` | Total requested CPUs per build/analysis, including the controller; null adds no limit |
| `total_limits.mem_gb` | Total requested memory in GB per build/analysis, including the controller; null adds no limit |
| `per_job_resources.<job>` | Individual job `cpus`, `mem_gb` (GB), and Slurm `time` |

The same `total_limits.jobs` value applies to sample arrays and Snakemake.
Arrays count running tasks; Snakemake counts queued and running workers.
With `jobs: null`, arrays omit the `%N` throttle and Snakemake uses `jobs: unlimited`.
Finite CPU/memory budgets still constrain scheduling. `array_size` and existing
batch dependencies remain in effect, as do the cluster's own resource and submission limits.
CPU and memory budgets limit requested allocations, not measured usage, and apply
independently to each build or analysis. They do not coordinate separate runs.

For each sample stage, concurrency is the minimum of `total_limits.jobs`,
`floor(total_limits.cpus / job.cpus)`, and
`floor(total_limits.mem_gb / job.mem_gb)`, omitting null limits. Thus a 512 GB
budget permits four 128 GB assembly tasks. Sample stages and the controller
remain sequential, using the existing Slurm dependencies.

Before scheduling Snakemake workers, the controller's CPU and memory requests
are subtracted from finite budgets. Queued workers also occupy budget until
completion. Effective rule requests, including per-job overrides, are charged
against the remaining budgets. A job that cannot fit fails with an error;
its CPU or memory request is never silently reduced to fit the budget.

`per_job_resources` combines sample stages (`assembly`, `busco`, `quant`), the
`controller`, and Snakemake rule names such as `odb_map`. For workers, explicit
config overrides take precedence over rule definitions. Unspecified worker
resources retain their rule definitions, with internal fallbacks of 8 GB memory
and one day; CPU counts retain each rule's thread count (one when unspecified).
`analysis.yaml` explicitly lists defaults for major worker steps so they can
be adjusted in one place. Entries for optional jobs do not enable those analyses.

`partition: null` uses the cluster default; use `sinfo` to check available names.
Account selection is left to the cluster; submission commands omit `--account`,
and the Snakemake profiles disable automatic account selection.
Keep job counts within the site's submission limit and `array_size` below
`MaxArraySize`. Slurm still determines when resources are available.

With `array_size: 1000`, `total_limits.jobs: 64`, and no tighter CPU/memory
budget, a full batch is equivalent to
`#SBATCH --array=1-1000%64`: 1000 sample tasks, at most 64 running at once.
Larger inputs are split into batches, with each batch waiting for the preceding
one to succeed. For example, 2500 pending samples in the assembly stage run as:

| Batch | Sample indices | Maximum running tasks | Starts after |
| --- | --- | --- | --- |
| 1 | 1–1000 | 64 | Resources become available |
| 2 | 1001–2000 | 64 | All tasks in batch 1 succeed |
| 3 | 2001–2500 | 64 | All tasks in batch 2 succeed |

All batches are submitted up front, linked with `--dependency=afterok:<job_id>`.
The chain keeps the sample-task concurrency limit across batches; a slow final
task holds back the next batch even when other slots are idle. A failed task
prevents dependent batches from starting. The next stage starts after all batches
of the preceding stage succeed. Retries submit only pending indices within each
batch. Splitting accommodates cluster array-size limits; the sequential batch
execution is a workflow scheduling choice.

Use `time` throughout `per_job_resources`.
It is the walltime limit of each job, using a quoted Slurm duration:
`time: "1-00:00:00"` means one day and `time: "02:00:00"` means two hours.
A rule's explicit `time` overrides the internal one-day fallback. Sample-stage
and controller limits are specified alongside rule limits in `per_job_resources`.
These durations are individual job limits, not a time budget for the entire workflow.

Prefer `"days-HH:MM:SS"` or `"HH:MM:SS"`. The other Slurm forms (`"minutes"`,
`"minutes:seconds"`, `"days-hours"`, and `"days-hours:minutes"`) are also accepted.
Always quote durations so YAML does not interpret colon-separated values as
numbers. Durations must be positive. Generated Snakemake profiles convert them
to the required internal `runtime` in minutes, rounding partial minutes up as
Slurm does. Replace old `runtime: 1440` settings with `time: "1-00:00:00"`.

`mem_gb` takes a positive integer and uses the same conversion as
`run_pipeline.sh --resources mem_gb=...`: multiply by 1000 for the internal
`mem_mb` value. For example, replace `mem_mb: 128000` with `mem_gb: 128`.
This preserves the previous memory request, including `--mem=128000M` for
sample-stage jobs. Generated and native Snakemake profiles still use `mem_mb`.

Resource requests are per job. For example, in `build.yaml`:

```yaml
slurm:
  total_limits:
    jobs: 64
    cpus: null
    mem_gb: 512
  per_job_resources:
    assembly: {cpus: 8, mem_gb: 256, time: "7-00:00:00"}
    odb_map: {cpus: 16, mem_gb: 192, time: "3-00:00:00"}
```

Apply resource edits to an existing preparation explicitly:

```bash
./run_build.sh submit --build results/angiosperm_leaf_20260925 --resources config/build.yaml
./run_analysis.sh submit --analysis results/angiosperm_leaf_20260925/downstream/c4_photosynthesis_20260929 --resources config/analysis.yaml
```

Only resources change; scientific settings remain frozen. Partial override files
can contain just the changed `total_limits` or `per_job_resources` entries.
Legacy `jobs`, `concurrency`, `stages`, `rules`, and `default_resources` settings
remain readable, including old retry files and snapshots. Old unequal job caps
are preserved until explicitly replaced by `total_limits.jobs`. Saved records
are not rewritten; the existing code/checksum verification still applies.
See [phylogeny](phylogeny.md#resources) and [dating](dating.md#resources) for their defaults.

## Re-running and recovery

Check `status`, `squeue`, and logs, then resubmit the same prepared build/analysis.
Completed work is reused; failed or missing work is retried. Active or unresolved
submissions block retries. A timed-out controller may leave workers running;
inspect those too.

New build IDs can reuse completed databases selected with `reuse_from`. New
analysis IDs share build products and reference caches, but do not automatically reuse
optional alignment, KO, or tree results from another analysis.
See [build failures](datasets.md#failure-and-recovery-behavior) and
[updating samples](datasets.md#updating-samples) for details.
