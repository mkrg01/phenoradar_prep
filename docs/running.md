# Running the workflow

[Documentation](index.md) · [Configuration](configuration.md)

## Installation and normal execution

Use Linux x86-64/Bash with the [requirements](../README.md#requirements).
Compute nodes need host Snakemake and `singularity` on `PATH`; Apptainer must
provide its `singularity` compatibility command.

The workflow image follows [VERSION](../VERSION). Use a published release and
prepare its image before the first mapping/analysis job or DAG dry-run:

```bash
./run_pipeline.sh --prepare-container --cores 1 --resources mem_gb=4
```

Repeat after changing releases; cached images are reused. The repository is
mounted automatically, so keep imported products inside it. Then follow the
[build/analysis commands](datasets.md) from metadata through completed outputs.

`plan` and `prepare` run locally. `submit` validates inputs, submits Slurm jobs,
and returns without waiting for completion; `submit --dry-run` writes/previews
submission scripts without submitting jobs. It is not a full Snakemake DAG dry-run.
Use `status`, `squeue`, and the prepared run's `jobs/logs/` to inspect progress.

Assembly, BUSCO, and quantification use species arrays with `afterok` dependencies.
Mapping and analysis each use a controller; Snakemake submits individual rules as
separate workers. Each array task, controller, and worker has its own time limit.
Status polling uses `squeue`, starting at 10 seconds; long-running jobs can increase
the interval automatically.

## Pilot run

For upstream inspection, add `--species-list pilot.txt` to build `submit`, then
submit without the list to finish the build. For a smaller analysis, set
`inputs.species_list` and `selection.species_list: true` in `analysis.yaml` before
preparing a new analysis. BUSCO filtering still applies to the candidate species.

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
| `timetree` | [treePL dating](dating.md) |
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

`partition: null` uses the cluster default. Check available names with `sinfo`
before choosing an explicit partition.

`jobs` is not a CPU count. Slurm starts only the jobs that fit the requested
CPU/memory allocations and site/account limits; the rest wait in the queue.
Set `jobs` within your site's submission limit; raising it does not reserve that many workers
at once. Assembly/BUSCO/quant arrays instead use `concurrency`.

`array_size` must be below the site's `MaxArraySize`. For example, increase assembly
memory/time or change ODB job sizing in `build.yaml`:

```yaml
slurm:
  stages:
    assembly: {cpus: 8, mem_mb: 256000, time: "7-00:00:00"}
  rules:
    odb_map: {cpus: 16, mem_mb: 192000, runtime: 4320}
```

Apply resource edits to an existing preparation explicitly:

```bash
./run_build.sh submit --build builds/build001 --resources config/build.yaml
./run_analysis.sh submit --analysis analyses/analysis001 --resources config/analysis.yaml
```

Only resources change; scientific settings remain frozen. Requests are per job,
not totals for all concurrent jobs. See [phylogeny](phylogeny.md#resources) and
[dating](dating.md#resources) for their defaults.

## Re-running and recovery

After checking the queue and logs, resubmit the same prepared build/analysis.
Completed work is reused within that run; failed or missing work is retried.
Active or unresolved submissions block overlapping retries. A timed-out controller
may leave workers active: inspect those before resubmitting.

New build IDs can reuse registered species products and ODB caches. New analysis
IDs share build products and reference caches, but do not automatically reuse
optional alignment, KO, or tree results from another analysis.
See [build failures](datasets.md#failure-and-recovery-behavior) for download and
assembly recovery, and [imports](datasets.md#updating-species-and-importing-existing-work)
for older results.

## Direct execution

A prepared analysis can run synchronously with a local resource budget:

```bash
./run_analysis.sh run --analysis analyses/analysis001 --local --cores 16 --mem-mb 192000
```

For low-level targets or a Snakemake DAG preview, use the **resolved** config,
not `config/build.yaml` or `config/analysis.yaml`:

```bash
./run_pipeline.sh --cores 16 --resources mem_gb=192 \
  --configfile analyses/analysis001/pipeline.yaml --dry-run -- all
```

Options precede `--`; targets follow it. `references`, `kegg_references`, and
`filter_species` are low-level targets, not `run_analysis.sh --target` values.
A resolved build config supports `mapping`; afterwards run
`run_build.sh complete --build builds/<id>` to validate and publish it.

The low-level launcher accepts `mem_mb` or whole decimal `mem_gb`. Inside a direct
Slurm allocation it derives the CPU/memory budget from Slurm and reserves 4 GB;
request one node, one task, and finite memory. Local runs should leave memory
outside the budget for the controller and other processes.

For host Conda execution, add `--software-deployment-method conda` to the low-level launcher.
