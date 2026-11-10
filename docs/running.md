# Installation and execution

[Documentation](index.md) · [Build and analysis guide](datasets.md)

## Installation and normal execution

Use Linux x86-64/Bash with Conda, Slurm or AGE, and Apptainer/Singularity. Compute nodes
need `snakemake` and `singularity` on `PATH`; Apptainer must provide its
`singularity` compatibility command.

From the repository root, create the host environment once and activate it in
each session:

```bash
conda env create -n phenoradar_prep -f environment.yaml
conda activate phenoradar_prep
```

For an existing installation, run `conda env update -n phenoradar_prep -f
environment.yaml` to add the AGE executor dependency.

The launchers use the active environment and do not activate it automatically.
Analysis tools run inside the workflow container. Use a published release and
prepare the image matching [VERSION](../VERSION):

```bash
./run_pipeline.sh --prepare-container --cores 1 --resources mem_gb=16
```

Prepare the image before the first mapping or analysis job, after changing
releases, or when the image cache has been cleared. GeneGalleon uses
[separate software](datasets.md#pinned-genegalleon-source-and-sif).

The first OCI-to-SIF conversion needs more memory than running an existing
image. Allocate at least 20 GB in Slurm to leave room for the launcher's 4 GB
overhead reservation. Compression uses at most the requested CPUs and one
quarter of the memory budget (capped at 4 GiB); temporary files stay on the
image-cache filesystem. Failed pulls cannot publish a partial image.
Job-specific limits require an unprivileged runtime installation; a setuid
installation must already have sufficiently strict administrator
`mksquashfs procs` and `mksquashfs mem` limits. System settings are never edited.

Continue with [building a database and running analyses](datasets.md).

The same launchers and configs work with both schedulers. `execution.scheduler:
auto` first checks the current allocation, then the available client commands
(`sbatch/squeue/sacct` or `qsub/qstat/qacct`). `run_pipeline.sh` also detects the
scheduler; explicit `--cores` runs locally. Inside an allocation it uses that
allocation unless called with `--cluster` by a build/analysis controller.

AGE discovers a single-node parallel environment and memory resource through
`qconf`. If detection is ambiguous, set `execution.scheduler: slurm` or `age`
(or `PHENORADAR_SCHEDULER`). Site resource overrides belong in the same config:

```yaml
execution:
  scheduler: auto
  age:
    parallel_environment: def_slot
    memory_resource: s_vmem
```

Those AGE overrides are normally unnecessary on HGC. Optional `age.queue` and
`age.extra_resources` select site resources. Job memory is a **total** request;
AGE converts it to per-slot memory. All slots must be on one node. AGE worker
completion requires `qacct`; disappearance from `qstat` does not imply success.
An ambiguous submission is recorded and blocks retry until reconciled, under
`jobs/` for controllers/arrays or `.snakemake/cluster_submissions/` for workers.

Slurm workers use greedy scheduling. A serialized `sbatch` gate limits actual
submissions to one per two seconds; the controller disables Snakemake's duplicate
rate-limit timers. Completion checks use `sacct`, so Slurm accounting must be
available and retain finished jobs. Temporary explicit rejections are retried up to three times. After an
ambiguous timeout, the launcher looks for its unique submission token in Slurm
before proceeding. If acceptance cannot be established, it stops and keeps the
receipt in `.snakemake/slurm_submissions/`; reconcile that receipt before retrying.

New downstream analyses create a compact, verified `build_inputs/` index.
Workers read only the selected product receipts rather than the full database.

## Resuming an unfinished v0.5.3 or v0.5.4 analysis with v0.5.5

After its controller and workers have ended, use the existing upgrade command:

```bash
python workflow/scripts/upgrade_analysis.py --analysis results/leaf/downstream/carnivory
python workflow/scripts/upgrade_analysis.py --analysis results/leaf/downstream/carnivory --apply
./run_analysis.sh submit --analysis results/leaf/downstream/carnivory
```

This narrowly audited migration backs up the saved code records and updates only
controller execution. It retains the saved tool image, all scientific settings,
outputs and output metadata; completed analyses cannot be migrated.
Read-only pending-work inspection uses the existing Conda/container cache paths
to recognize reusable outputs without installing or pulling environments.

## Resuming v0.5.7 after a completion-validation failure

Final native-output validation runs in the host environment, including `ete4`
for tree checks. With v0.5.8, update that environment first:

```bash
conda env update -n phenoradar_prep -f environment.yaml
conda activate phenoradar_prep
python workflow/scripts/upgrade_analysis.py --analysis results/leaf/downstream/carnivory
python workflow/scripts/upgrade_analysis.py --analysis results/leaf/downstream/carnivory --apply
./run_analysis.sh submit --analysis results/leaf/downstream/carnivory --dry-run
./run_analysis.sh submit --analysis results/leaf/downstream/carnivory
```

Wait for the previous controller and workers to end before migrating. The
upgrade saves the old code records and retains scientific settings, the saved
tool image, outputs and output metadata. Check that the preview reuses completed
workflow steps; whole-analysis completion still requires successful validation
and a new `completed.json`.

## Upgrading an unfinished v0.4.14 through v0.4.18 analysis

With a clean v0.4.19 checkout, audit compatible reuse before applying it.
For an existing analysis, explicitly opt in to the new 50% OG occupancy condition:

```bash
python workflow/scripts/upgrade_analysis.py --analysis results/leaf/downstream/carnivory --min-sample-fraction 0.5
python workflow/scripts/upgrade_analysis.py --analysis results/leaf/downstream/carnivory --min-sample-fraction 0.5 --apply
```

Run this verification with an allocated controller-sized CPU/memory budget.
It requires an idle repository, checks the pinned predecessor and unchanged
scientific code/environments, verifies imported artifacts (including decoded
content when legacy plain products were compressed), and backs up the old
records under the next `jobs/code_migration_*/` directory. Omitting
`--min-sample-fraction` preserves the saved condition; legacy analyses retain all
OGs with an explicit `0`. Condition changes update the frozen settings and are
recorded with before/after values. Only verified import-job metadata is reset.
A pinned compatibility launcher retains input/mtime/parameter/code rerun
checks while allowing the verified unchanged tool environments to be reused.
Use the normal submission command afterward; retain four FAMSA threads
per OG to reuse existing alignments. New batch resource overrides use
`per_job_resources.align_orthogroup_batch` (12 CPUs/32 GB by default).

## Targets

To run a branch in an existing analysis:

```bash
./run_analysis.sh submit --analysis results/leaf/downstream/carnivory --target phylogeny
```

Targets include missing prerequisites within the analysis. They use the saved
settings and never rebuild upstream sample products.

| Target | Work requested |
| --- | --- |
| `all` (default) | Merge OG/per-sequence expression, run enabled branches, validate native outputs, and record completion |
| `alignments` | [All-copy OG alignments](alignments.md) and per-sequence TPM |
| `kegg` | [KO annotation and expression](kegg.md) |
| `phylogeny_prepare` | BUSCO input audit, outgroup, and marker plan |
| `phylogeny` | Infer trees selected by `phylogeny.trees` |
| `phylogeny_calibrations` | Trees and [TimeTree calibrations](dating.md#timetree-calibrations), without dating |
| `timetree` | [LSD2 dating](dating.md) |
| `taxonomy_check` | [MonoPhy review](taxonomy_check.md) |
| `contrast_pairs` | [Trait pairs](contrast_pairs.md) on selected trees |

Dating, taxonomy checks, and contrast pairs require their enabled flag even
when targeted explicitly. `alignments` and `kegg` work without enabling their
inclusion in `all`. See [build execution](datasets.md#build-a-reusable-database)
for the complete per-sample native workflow and database publication.

## Resource budgets

Edit `execution` in the relevant build or analysis config:

| Setting | Controls |
| --- | --- |
| `scheduler` | `auto` (default), `slurm`, or `age` |
| `partition` | Slurm partition; `null` uses the cluster default |
| `total_limits.jobs` | Running array tasks or queued/running Snakemake workers; excludes the controller |
| `total_limits.cpus`, `total_limits.mem_gb` | Total requested CPUs/memory for one build or analysis, including its controller |
| `per_job_resources.<job>` | Per-job `cpus`, `mem_gb`, and `time` |

A `null` total limit adds no workflow cap. Limits apply independently to each
run and measure requested allocations, not actual usage. Cluster limits still
apply.

The build configs, including angiosperm, default to `total_limits.jobs: null`.
New analyses keep
finite CPU/memory budgets (64 CPUs/256 GB by default), so worker concurrency
follows each job's allocation rather than a fixed worker count. Prepared analyses
keep their saved settings; use resource overrides to change a later submission.

```yaml
execution:
  total_limits:
    jobs: 64
    cpus: null
    mem_gb: 512
  per_job_resources:
    sample: {cpus: 8, mem_gb: 128, time: "16-00:00:00"}
    odb_map: {cpus: 16, mem_gb: 128, time: "2-00:00:00"}
```

Here the memory budget allows at most four sample tasks at once. For Snakemake
workers, the controller's allocation is subtracted first; queued workers also
occupy budget. A job that cannot fit fails with an error rather than receiving
fewer resources.

`per_job_resources.sample` sizes the entire assembly/BUSCO/quantification job.
The supplied config requests 4 CPUs, 128 GB, and 16 days per sample; the allocation
remains reserved throughout processing.

Samples run as array tasks, split into sequential batches to fit the scheduler's
array limits and the workflow budget. Preparation reads these limits through
`scontrol show config` (Slurm) or `qconf -sconf` (AGE) before
submitting arrays. Failed samples do not block other samples; database publication
waits for all of them to succeed.

`per_job_resources.controller` sizes preparation and the controllers
that manage mapping and analysis jobs. The default build controller requests
four CPUs. Checksum verification uses up to four threads, limited to allocated
CPUs.

OG expression aggregation processes up to 50 runs per job. Set
`per_job_resources.aggregate_tpm` for each batch (default: one CPU and 2 GB).
Logs are named `orthogroups/expression/batch_*.log`; retries reuse completed runs.

`merge_tpm` combines the run tables in one job and reads each table once. It
retains compact expression columns and writes wide tables one run at a time.
The default request is one CPU and 32 GB; use `per_job_resources.merge_tpm` to
override it.

Resource overrides take precedence over rule defaults. Enable optional
branches in the analysis settings.

Memory uses positive integer GB (`mem_gb: 128` becomes `128000` MB). Both backends
accept the same quoted
walltimes, such as `"1-00:00:00"` for one day or `"02:00:00"` for two hours.
Time limits apply to each job, including the controller.

To change resources for an existing run:

```bash
./run_build.sh submit --build results/leaf --resources config/build.yaml
./run_analysis.sh submit --analysis results/leaf/downstream/carnivory --resources config/analysis.yaml
```

Only the supplied `execution` settings change; scientific settings stay frozen.
Existing configs and saved builds using `slurm` remain supported; do not specify
both `execution` and `slurm` in one file.
Partial resource files may contain just the changed limits or job entries.

GeneGalleon derives its tool limits from the actual allocation. Worker logs
show the allocation and tool budgets; retry resource overrides also change the
assembler's limits.

## Re-running and recovery

Use the run's `status` command, `squeue`/`qstat`, and `jobs/logs/` to inspect progress.
See [failure recovery](datasets.md#failure-and-recovery-behavior) for retry commands
and [database reuse](datasets.md#reusing-completed-databases) for new builds.
