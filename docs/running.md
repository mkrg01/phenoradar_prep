# Installation and execution

[Documentation](index.md) · [Build and analysis guide](datasets.md)

## Installation and normal execution

Use Linux x86-64/Bash with Conda, Slurm, and Apptainer/Singularity. Compute nodes
need `snakemake` and `singularity` on `PATH`; Apptainer must provide its
`singularity` compatibility command.

From the repository root, create the host environment once and activate it in
each session:

```bash
conda env create -n phenoradar_prep -f environment.yaml
conda activate phenoradar_prep
```

The launchers use the active environment and do not activate it automatically.
Analysis tools run inside the workflow container. Use a published release and
prepare the image matching [VERSION](../VERSION):

```bash
./run_pipeline.sh --prepare-container --cores 1 --resources mem_gb=4
```

Prepare the image before the first mapping or analysis job, after changing
releases, or when the image cache has been cleared. GeneGalleon uses
[separate software](datasets.md#pinned-genegalleon-source-and-sif).

Continue with [building a database and running analyses](datasets.md).

## Targets

To run a branch in an existing analysis:

```bash
./run_analysis.sh submit --analysis results/leaf/downstream/carnivory --target phylogeny
```

Targets include missing prerequisites within the analysis. They use the saved
settings and never rebuild upstream sample products.

| Target | Work requested |
| --- | --- |
| `all` (default) | Merge OG expression/QC, run enabled branches, collect PhenoRadar inputs |
| `alignments` | [All-copy OG alignments](alignments.md) |
| `kegg` | [KO annotation and expression](kegg.md) |
| `phylogeny_prepare` | BUSCO input audit, outgroup, and marker plan |
| `phylogeny` | Infer trees selected by `phylogeny.trees` |
| `phylogeny_calibrations` | Trees and [TimeTree calibrations](dating.md#timetree-calibrations), without dating |
| `timetree` | [LSD2 dating](dating.md) |
| `taxonomy_check` | [MonoPhy review](taxonomy_check.md) |
| `contrast_pairs` | [Trait pairs](contrast_pairs.md) on selected trees |
| `phenoradar_inputs` | [Collect completed results](phenoradar_inputs.md), without starting analyses |

Dating, taxonomy checks, and contrast pairs require their enabled flag even
when targeted explicitly. `alignments` and `kegg` work without enabling their
inclusion in `all`. See [build execution](datasets.md#build-a-reusable-database)
for the complete per-sample native workflow and database publication.

## Resource budgets

Edit `slurm` in the relevant build or analysis config:

| Setting | Controls |
| --- | --- |
| `partition` | Slurm partition; `null` uses the cluster default |
| `total_limits.jobs` | Running array tasks or queued/running Snakemake workers; excludes the controller |
| `total_limits.cpus`, `total_limits.mem_gb` | Total requested CPUs/memory for one build or analysis, including its controller |
| `per_job_resources.<job>` | Per-job `cpus`, `mem_gb`, and `time` |

A `null` total limit adds no workflow cap. Limits apply independently to each
run and measure requested allocations, not actual usage. Cluster limits still
apply.

```yaml
slurm:
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

Samples run as array tasks, split into sequential batches to fit
[Slurm's array limits](https://slurm.schedmd.com/job_array.html) and the workflow
budget. Preparation must read these limits through `scontrol show config` before
submitting arrays. Failed samples do not block other samples; database publication
waits for all of them to succeed.

`per_job_resources.controller` sizes Slurm preparation and the controllers
that manage mapping and analysis jobs. The default build controller requests
four CPUs. Checksum verification uses up to four threads, limited to allocated
CPUs.

OG expression aggregation processes up to 50 runs per Slurm job. Set
`per_job_resources.aggregate_tpm` for each batch (default: one CPU and 2 GB).
Logs are named `orthogroups/expression/batch_*.log`; retries reuse completed runs.

`merge_tpm` combines the run tables in one job and reads each table once. It
retains compact expression columns and writes wide tables one run at a time.
The default request is one CPU and 32 GB; use `per_job_resources.merge_tpm` to
override it.

Resource overrides take precedence over rule defaults. Enable optional
branches in the analysis settings.

Memory uses positive integer GB (`mem_gb: 128` becomes `128000` MB). Quote Slurm
walltimes, such as `"1-00:00:00"` for one day or `"02:00:00"` for two hours.
Time limits apply to each job, including the controller.

To change resources for an existing run:

```bash
./run_build.sh submit --build results/leaf --resources config/build.yaml
./run_analysis.sh submit --analysis results/leaf/downstream/carnivory --resources config/analysis.yaml
```

Only the supplied `slurm` settings change; scientific settings stay frozen.
Partial resource files may contain just the changed limits or job entries.

GeneGalleon derives its tool limits from the actual Slurm allocation. Worker logs
show the allocation and tool budgets; retry resource overrides also change the
assembler's limits.

## Re-running and recovery

Use the run's `status` command, `squeue`, and `jobs/logs/` to inspect progress.
See [failure recovery](datasets.md#failure-and-recovery-behavior) for retry commands
and [database reuse](datasets.md#reusing-completed-databases) for new builds.
