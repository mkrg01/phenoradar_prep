# Running the workflow

[Documentation](index.md) · [Configuration](configuration.md)

## Installation and normal execution

Use Slurm on Linux x86-64/Bash with the [requirements](../README.md#requirements).
Compute nodes need host Snakemake and `singularity` on `PATH`; Apptainer must
provide its `singularity` compatibility command.

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

`plan` and `prepare` run locally; `submit` validates inputs, submits Slurm jobs,
and returns without waiting for completion. `submit --dry-run` previews submission
scripts without submitting jobs; it does not check the full Snakemake DAG.
Use `status`, `squeue`, and the prepared run's `jobs/logs/` to inspect progress.

Assembly, BUSCO, and quantification use sample arrays with `afterok` dependencies.
Mapping and analysis use a controller with separate Snakemake worker jobs.
Each array task, controller, and worker has its own resource and time limits.

## Targets

Use `run_analysis.sh submit --analysis analyses/<analysis> --target <target>`.
Targets include missing prerequisites within the analysis; they never rebuild
upstream species products. Edit branch settings before `prepare`.

| Analysis target | Work requested |
| --- | --- |
| `all` (default) | OG expression/QC and enabled branches, followed by PhenoRadar collection |
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
Build endpoints are documented [separately](datasets.md#build-through-mapping).

## Resource budgets

Edit the relevant config's `slurm` section:

| Setting | Controls |
| --- | --- |
| `partition`, `account` | Site allocation |
| `concurrency`, `array_size` (build only) | Concurrent species tasks and array batch size |
| `stages.assembly`, `.busco`, `.quant` (build only) | Per-species `cpus`, `mem_mb`, and Slurm `time` |
| `stages.controller` | Controller `cpus`, `mem_mb`, and `time` |
| `jobs` | Maximum outstanding Snakemake worker jobs (default 64), including queued jobs |
| `rules.<rule>` | Per-rule `cpus`, `mem_mb`, and `runtime` in minutes |

`partition: null` uses the cluster default; use `sinfo` to check available names.
Keep `jobs` within the site's submission limit and `array_size` below `MaxArraySize`.
`jobs` counts queued/running workers, not CPUs; arrays use `concurrency` instead.
Resource requests are per job. For example, in `build.yaml`:

```yaml
slurm:
  stages:
    assembly: {cpus: 8, mem_mb: 256000, time: "7-00:00:00"}
  rules:
    odb_map: {cpus: 16, mem_mb: 192000, runtime: 4320}
```

Apply resource edits to an existing preparation explicitly:

```bash
./run_build.sh submit --build builds/angiosperm_leaf_20260925 --resources config/build.yaml
./run_analysis.sh submit --analysis analyses/analysis001 --resources config/analysis.yaml
```

Only resources change; scientific settings remain frozen. See
[phylogeny](phylogeny.md#resources) and [dating](dating.md#resources) for their defaults.

## Re-running and recovery

Check `status`, `squeue`, and logs, then resubmit the same prepared build/analysis.
Completed work is reused; failed or missing work is retried. Active or unresolved
submissions block retries. A timed-out controller may leave workers running;
inspect those too.

New build IDs can reuse registered species products and ODB caches. New analysis
IDs share build products and reference caches, but do not automatically reuse
optional alignment, KO, or tree results from another analysis.
See [build failures](datasets.md#failure-and-recovery-behavior) and
[updating samples](datasets.md#updating-samples) for details.
