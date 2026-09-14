# Running the workflow

[Documentation](index.md) · [Configuration](configuration.md)

Run commands from the repository root after configuring the [inputs](inputs.md).
The launcher runs Snakemake locally, either on the current host or inside one
Slurm allocation. It prints shell commands and reruns incomplete jobs automatically.

## Installation

Use Linux with Bash and Conda:

```bash
conda env create -n phenoradar-workflow -f environment.yaml
conda activate phenoradar-workflow
```

This installs Snakemake. Processing tools have separate environments in
`workflow/envs/`; direct execution needs `--software-deployment-method conda`
to create and use them. The launcher enables Conda automatically inside Slurm.
Definitions pin the main packages but are not complete dependency lockfiles.

If using an existing software installation, expose the workflow's fixed commands
on `PATH` and omit the deployment option during direct execution.
`SNAKEMAKE_BIN` can select the Snakemake executable. The launcher stores its
Snakemake cache in `.cache/` through `XDG_CACHE_HOME`.

Phylogeny prepares the [ASTRAL tool](phylogeny.md#setup-and-execution) on first use.
For container releases and Apptainer execution, see [containers](containers.md).
Missing database snapshots are [prepared automatically](references.md) when a
requested branch needs them. First use therefore needs download access unless
references and software have been prepared locally.

## Targets

Put targets after `--`, with every option before it:

```bash
./run_pipeline.sh --software-deployment-method conda \
  --configfile config/mydata.yaml --cores 2 --resources mem_gb=16 -- prepare
```

Explicit analysis targets work even when their `enabled` flag is false.
Prerequisites are scheduled automatically, including missing or outdated tree
inference for `taxonomy_audit` and unfiltered `phylogeny_contrast_pairs`.

| Target | Work requested |
| --- | --- |
| `all` (default) | OG expression/QC and branches enabled by `alignment`, `kegg`, `phylogeny`, and `contrast` |
| `prepare` | Metadata selection, PhenoRadar species metadata, run record, sample/chunk manifests |
| `references` | OrthoDB snapshot only |
| `kegg_references` | KOfam/KEGG snapshot only; no assemblies required |
| `proteins` | CDS translation for selected species |
| `mapping` | ODB mapping and merged gene-to-OG index |
| `alignments` | [All-copy OG alignments](alignments.md), including mapping |
| `kegg` | [KO annotation and original-TPM sums](kegg.md), independent of ODB |
| `phylogeny_prepare` | BUSCO input audit, outgroup resolution, and marker plan |
| `phylogeny` | [BUSCO species-tree inference](phylogeny.md); also dating/audit when their flags are enabled |
| `phylogeny_calibrations` | Species-tree inference and [TimeTree calibration retrieval](dating.md#timetree-calibrations), without dating |
| `timetree` | Species-tree inference and [LSD2 dating](dating.md) using the selected calibration source |
| `taxonomy_audit` | [MonoPhy review](taxonomy_audit.md) of full/phenotyped species trees |
| `contrast_pairs` | [Representative selection, inference, and trait pairs](contrast_pairs.md#representative-analysis) |
| `phylogeny_contrast_pairs` | [Trait pairs from full/phenotyped trees](contrast_pairs.md#pairs-from-full-or-phenotyped-trees); with exclusions, requires completed results and uses the filtered export |
| `phenoradar_metadata` | Minimal species metadata, useful for backfilling older results |
| `filter_species` | [Export completed results after exclusions](species_filter.md); does not start producer analyses |
| `phenoradar_inputs` | [Validate and collect completed inputs](phenoradar_inputs.md); does not start producer analyses |

The full/phenotyped phylogeny targets follow `phylogeny.species_sets`.
`contrast.enabled` adds only the separate representative analysis to `all`.
Filtering and PhenoRadar collection are always manual targets.

## Slurm

Activate the workflow environment and create `logs/` before submitting. Replace
`YOUR_PARTITION` and add your cluster's account option if needed:

```bash
mkdir -p logs
sbatch --partition=YOUR_PARTITION \
  run_pipeline.sh --configfile config/mydata.yaml
```

The script requests one node, one task, 16 CPUs, 192 GiB, and up to 21 days by
default. All steps share that allocation. Set allocation options before the
script name and workflow options after it. For example, preparation can use:

```bash
sbatch --partition=YOUR_PARTITION --cpus-per-task=2 --mem=16G --time=01:00:00 \
  run_pipeline.sh --configfile config/mydata.yaml -- prepare
```

The first taxonomy build may need a longer allocation. With a prepared snapshot,
8 GiB fits the declared preparation budgets. Avoid concurrent submissions writing
the same outputs.

Monitor or stop a job using its printed ID:

```bash
squeue -u "$USER"
tail -f logs/pipeline-JOB_ID.log
scancel JOB_ID
```

Resubmit the same command to resume after cancellation or interruption. The
script requests an interrupt shortly before its time limit and preserves
Snakemake's exit status as the batch job's exit status.

## Direct execution

Supply a CPU and memory budget for the current host:

```bash
./run_pipeline.sh --software-deployment-method conda \
  --configfile config/mydata.yaml --cores 16 --resources mem_gb=192
```

Leave physical memory for Snakemake and other processes beyond the scheduling
budget. The `#SBATCH` lines are comments in this mode and reserve no resources.
Inside an existing Slurm allocation (`SLURM_JOB_ID` is set), the launcher instead
uses that allocation's resources and enables Conda. Run from the repository root
in that case.

## Resource budgets

Rule threads and memory settings describe one job. The launcher budget limits
the sum of concurrent jobs. For example, two default ODB chunks need 32 CPUs
and 384 GB of workflow memory; one needs 16 CPUs and 192 GB. Increase the
allocation to permit more concurrency:

```bash
sbatch --partition=YOUR_PARTITION --cpus-per-task=32 --mem=384G \
  run_pipeline.sh --configfile config/mydata.yaml
```

Rule thread counts are capped by available cores. Rule memory values are
scheduling estimates, while Slurm enforces total allocation memory. Defaults
should be checked against pilot benchmarks. An allocation must fit the largest
step and last long enough for all chunks; its resources remain reserved during
lighter stages too.

Configuration `mem_gb` values and direct `--resources mem_gb=...` use positive
whole decimal GB (1 GB = 1000 MB), converted internally to Snakemake `mem_mb`.
Slurm `--mem=192G` uses GiB. The launcher converts Slurm's units and leaves 4 GB
for Snakemake and environment management; that default allocation provides
about 202.2 GB for jobs, enough for one 192 GB ODB chunk.

The launcher enforces the local executor and, under Slurm, allocation-derived
budgets. Request one node, one task, and finite memory with `--mem` or
`--mem-per-cpu`; set total resources when requesting the allocation.

## Prepare and inspect

Run `prepare` to select species and create manifests before ODB reference
preparation or mapping. Selection is a checkpoint, so a dry-run before preparation
may not show every downstream job. Inspect the plan again after preparation:

```bash
./run_pipeline.sh --software-deployment-method conda \
  --configfile config/mydata.yaml --cores 16 --resources mem_gb=192 --dry-run
```

Review `metadata/selection.json`, `samples.tsv`, and `busco_completeness.svg`
under your analysis directory before scaling up.

## Pilot run

After preparation, choose a few eligible IDs from
`results/<analysis>/metadata/species_high_busco.txt`. For the default `full`
analysis, this gives a starting list to review:

```bash
head -n 3 results/full/metadata/species_high_busco.txt > input/pilot_species.txt
sbatch --partition=YOUR_PARTITION --cpus-per-task=8 --mem=80G \
  run_pipeline.sh --configfile config/mydata.yaml config/pilot.yaml
```

`config/pilot.yaml` uses two-species ODB chunks, retains work directories, and
writes to `results/pilot/`. The species list is user-supplied and ignored by Git.
Check mapping quality, runtime, disk use, and peak memory before the full run.

## Re-running and recovery

Snakemake tracks declared inputs, parameters, code, and environment definitions.
Changing only an abundance file recalculates expression for that run and the
combined tables. Changing selected species updates the relevant manifests and
analyses. Final merged tables use only current species/runs, even when older
per-run files remain on disk.

ODB jobs retain interrupted work under
`work/<analysis>/orthogroups/mapping/<chunk>/<fingerprint>/`. The fingerprint
includes protein contents, reference, mapping options, software records, and
worker code. Matching work can resume; changed inputs use a new directory.
Chunk membership, `odb.threads`, and `odb.batch_size` affect this identity;
`odb.mem_gb` alone does not.

Successful ODB jobs copy native results and logs out before removing their work
directory, unless `odb.keep_work: true`. Work for other fingerprints is retained.
Replacing installed ODB software in place may require an explicit rerun:

```bash
sbatch --partition=YOUR_PARTITION run_pipeline.sh \
  --configfile config/mydata.yaml --forcerun odb_map
```

Reference updates have separate [snapshot refresh procedures](references.md).
Use the [migration guide](migration.md) when resuming results from older layouts
or configurations.
