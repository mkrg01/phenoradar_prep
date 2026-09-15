# Running the workflow

[Documentation](index.md) · [Configuration](configuration.md)

## Installation

Install the [requirements](../README.md#requirements) on Linux with Bash and make
`snakemake` and `singularity` available on `PATH`, including compute nodes.
`SNAKEMAKE_BIN` can select another Snakemake executable. No dedicated host Conda
environment is required.

Use a published release checkout, edit [config/config.yaml](../config/config.yaml)
for your [inputs](inputs.md), and leave `container_image: auto` for automatic
image selection. The examples below use this configuration file automatically.
See [container setup](containers.md)
for version selection, external bind mounts, or native execution. First use needs network
access for the image and missing [references](references.md).

## Slurm

Adjust the `#SBATCH` lines in `run_pipeline.sh` for your cluster, especially
partition, account, CPUs, memory, and time. Submit from the repository root:

```bash
sbatch run_pipeline.sh
```

The workflow automatically uses the allocated CPUs and memory, so there is no
need to pass `--cores` or `--resources`. Omit the target to run `all`, including
preparation and enabled analyses.

All jobs run inside one node/task allocation. The script defaults to 16 CPUs,
192 GiB, and 21 days on `debug`. Override these with `sbatch` options before the
script name; workflow options go after it:

```bash
sbatch --cpus-per-task=2 --mem=16G --time=01:00:00 \
  run_pipeline.sh -- prepare
```

Slurm writes standard output to `pipeline-<job_id>.out` and standard error to
`pipeline-<job_id>.err` in the repository root. Per-step logs are written under
`logs/<run_name>/`, which Snakemake creates automatically.

Monitor with `squeue -u "$USER"` and the job's `.out` and `.err` files; stop with
`scancel JOB_ID`. Resubmit the same command to resume. Avoid concurrent jobs
writing the same results.

## Direct execution

Run from the repository root with a CPU and memory budget:

```bash
./run_pipeline.sh --cores 16 --resources mem_gb=192
```

Put all options before `--` and targets after it. Omit the target to run `all`.
For a separate configuration file, see [configuration](configuration.md#loading-settings-and-paths).

## Prepare and inspect

To inspect species selection and manifests before the full analysis, run `prepare`:

```bash
./run_pipeline.sh --cores 2 --resources mem_gb=16 -- prepare
```

Review `selection.json`, `samples.tsv`, and `busco_completeness.svg` in
`results/<run_name>/metadata/`. Then inspect downstream work with `--dry-run`.
Selection is a checkpoint, so a dry-run before preparation can be incomplete.

## Pilot run

Choose a few eligible species, then use the supplied pilot override:

```bash
head -n 3 results/run001/metadata/species_high_busco.txt > input/pilot_species.txt
sbatch --cpus-per-task=8 --mem=80G \
  run_pipeline.sh --configfile config/pilot.yaml \
  --set-threads odb_map=8 --set-resources odb_map:mem_mb=64000
```

[config/pilot.yaml](../config/pilot.yaml) overrides the species list and run name
in `config/config.yaml`: it uses `input/pilot_species.txt` and writes
to `results/pilot/`. Other settings come from `config/config.yaml`.
Review the species list and check mapping
quality, runtime, disk use, and peak memory before a full run.

## Targets

Explicit analysis targets work even when their `enabled` flag is false and
schedule missing prerequisites automatically.

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
| `phylogeny` | [BUSCO species-tree inference](phylogeny.md); also dating/taxonomy checks when their flags are enabled |
| `phylogeny_calibrations` | Species-tree inference and [TimeTree calibration retrieval](dating.md#timetree-calibrations), without dating |
| `timetree` | Species-tree inference and [LSD2 dating](dating.md) using the selected calibration source |
| `taxonomy_check` | [MonoPhy review](taxonomy_check.md) of full/phenotyped species trees |
| `contrast_pairs` | [Representative selection, inference, and trait pairs](contrast_pairs.md#representative-analysis) |
| `phylogeny_contrast_pairs` | [Trait pairs from full/phenotyped trees](contrast_pairs.md#pairs-from-full-or-phenotyped-trees); with exclusions, requires completed results and uses the filtered export |
| `phenoradar_metadata` | Minimal species metadata, useful for backfilling older results |
| `filter_species` | [Export completed results after exclusions](species_filter.md); does not start producer analyses |
| `phenoradar_inputs` | [Automatically collect available completed results](phenoradar_inputs.md); no collection settings or producer analyses |

The full/phenotyped phylogeny targets follow `phylogeny.species_sets`.
`contrast.enabled` adds only the separate representative analysis to `all`.
Filtering and PhenoRadar collection are always manual targets.

## Resource budgets

CPU and memory defaults live in the rules. Workflow configuration files contain
no CPU or memory settings. Rule resources apply to one job; the launcher budget
limits concurrent jobs.
A default ODB chunk requests 16 CPUs and 192 GB. Two concurrent chunks therefore
need 32 CPUs and 384 GB. Fit the largest step and check estimates with a pilot.

| Rule | Job unit | Default threads | Default memory (GB) |
| --- | --- | --- | --- |
| `odb_map` | Mapping chunk | 16 | 192 |
| `align_orthogroup` | OG | 4 | 8 |
| `annotate_kofam` | Species | 4 | 8 |
| `check_taxonomy` | Species set | 1 | 8 |

See [phylogeny resources](phylogeny.md#resources) for tree-inference defaults.
For an individual rule, use Snakemake's standard overrides:

```bash
./run_pipeline.sh --cores 16 --resources mem_gb=192 \
  --set-threads odb_map=8 --set-resources odb_map:mem_mb=64000 -- mapping
```

This requests 8 threads and 64 GB per ODB chunk within a total budget of 16 CPUs
and 192 GB. Chunks contain up to 100 species, defined by `make_manifests` in
`workflow/rules/odb.smk`. Internal batch size is four times the actual ODB thread
count (32 in this example), including any CPU cap applied by Snakemake.
Rule overrides use `mem_mb` (64000 MB = 64 GB). The launcher's total-budget option
`--resources mem_gb=...` uses positive whole decimal GB.
Slurm `--mem` uses GiB; the launcher converts units and reserves 4 GB for overhead.
Inside a Slurm allocation, its CPU/memory limits override direct launcher budgets.
Request one node, one task, and finite memory. For direct execution, leave
physical memory for Snakemake and other processes beyond the scheduling budget.

## Re-running and recovery

Rerun the same command after interruption. Snakemake reuses completed jobs when
their inputs, settings, and code are unchanged;
abundance-only changes recalculate expression without remapping proteins.

ODB work remains under `work/<run_name>/orthogroups/mapping/` after success or
failure. Matching inputs, settings, software, and code reuse the same directory,
allowing ODB to resume its internal steps. Changed inputs use separate work.
Native results and logs are also published under `results/<run_name>/`.
After replacing ODB software in place, explicitly rerun it with `--forcerun odb_map`.
Use a new `run_name` for fresh work.

KofamScan also retains work files under `work/<run_name>/kegg/`. Valid completed
species annotations are reused; failed annotations restart while retaining
their previous attempts. Other completed workflow steps follow Snakemake's
usual reuse rules.

See [reference updates](references.md) and [migration](migration.md) when changing
snapshots or reusing older results.
