# Running the workflow

[Documentation](index.md) · [Configuration](configuration.md)

## Installation

Use a published release checkout with the [requirements](../README.md#requirements)
on Linux/Bash, and edit [config/config.yaml](../config/config.yaml).
See [container setup](containers.md) for compute-node tools, images, and bind mounts.
First use needs network for the image and missing [references](references.md).
`SNAKEMAKE_BIN` can select another Snakemake executable.

## Slurm

Adjust the `#SBATCH` lines in `run_pipeline.sh` for your cluster (partition,
account, CPUs, memory, time), then submit from the repository root:

```bash
sbatch run_pipeline.sh
```

Jobs share one node/task allocation and automatically use its CPUs/memory.
Override Slurm settings before the script name:

```bash
sbatch --cpus-per-task=32 --mem=256G --time=7-00:00:00 run_pipeline.sh
```

Workflow options go after the script name. Monitor `pipeline-<job_id>.out` / `.err`
in the repository root and per-step logs in `logs/<run_name>/`.
Avoid concurrent jobs writing the same results.

## Direct execution

From the repository root, supply a CPU and memory budget:

```bash
./run_pipeline.sh --cores 16 --resources mem_gb=192
```

Put options before `--` and targets after it. Omitting the target runs `all`.
`config/config.yaml` loads automatically; use `--configfile` for
[overrides](configuration.md#loading-settings-and-paths).

## Pilot run

List a few candidate species IDs in `input/pilot_species.txt`, one per line.
[BUSCO filtering](inputs.md#species-selection) still applies. Then run:

```bash
sbatch --cpus-per-task=8 --mem=80G \
  run_pipeline.sh --configfile config/pilot.yaml \
  --set-threads odb_map=8 --set-resources odb_map:mem_mb=64000
```

[config/pilot.yaml](../config/pilot.yaml) overrides only species selection and
run name. Check selection/mapping QC in `results/pilot/`, runtime, disk use,
and peak memory before a full run.

## Targets

Explicit analysis targets work even when their `enabled` flag is false and
schedule missing prerequisites automatically.

| Target | Work requested |
| --- | --- |
| `all` (default) | OG expression/QC and enabled alignment, KEGG, phylogeny, and contrast branches |
| `references`, `kegg_references` | [OrthoDB or KOfam/KEGG snapshots](references.md) |
| `proteins` | CDS translation |
| `mapping` | ODB mapping and merged gene-to-OG index |
| `alignments` | [All-copy OG alignments](alignments.md) |
| `kegg` | [KO annotation and original-TPM sums](kegg.md) |
| `phylogeny_prepare` | BUSCO input audit, outgroup, and marker plan |
| `phylogeny` | [Species trees](phylogeny.md), plus enabled dating/taxonomy checks |
| `phylogeny_calibrations` | Trees and [TimeTree calibrations](dating.md#timetree-calibrations), without dating |
| `timetree` | Trees and [LSD2 dating](dating.md) |
| `taxonomy_check` | [MonoPhy review](taxonomy_check.md) |
| `contrast_pairs` | [Representative selection, inference, and pairs](contrast_pairs.md#representative-analysis) |
| `phylogeny_contrast_pairs` | [Pairs from full/phenotyped trees](contrast_pairs.md#pairs-from-full-or-phenotyped-trees); with exclusions, uses completed filtered results |
| `filter_species` | [Export completed results after exclusions](species_filter.md) |
| `phenoradar_inputs` | [Collect available completed results](phenoradar_inputs.md) |

Full/phenotyped targets follow `phylogeny.species_sets`. Filtering and collection
are manual targets that never start producer analyses.

## Resource budgets

Rule resources apply to one job; the launcher budget limits concurrent jobs.
Defaults live in the rules, not configuration files.

| Rule | Job unit | Threads | Memory (GB) |
| --- | --- | --- | --- |
| `odb_map` | Mapping chunk (up to 100 species) | 16 | 192 |
| `align_orthogroup` | OG | 4 | 8 |
| `annotate_kofam` | Species | 4 | 8 |
| `check_taxonomy` | Species set | 1 | 8 |

See [phylogeny resources](phylogeny.md#resources) for tree-inference defaults.
Two default ODB jobs need 32 CPUs and 384 GB. To reduce per-job requests:

```bash
./run_pipeline.sh --cores 16 --resources mem_gb=192 \
  --set-threads odb_map=8 --set-resources odb_map:mem_mb=64000 -- mapping
```

This allocates 8 threads and 64 GB per chunk within a 16-CPU/192-GB budget.

Rule overrides use `mem_mb`; the launcher uses positive whole decimal
`mem_gb`. Slurm `--mem` uses GiB; the launcher converts it and reserves 4 GB.
Inside Slurm, allocation limits override direct budgets; request one node, one
task, and finite memory. For direct runs, leave memory outside the budget for
Snakemake and other processes.

## Re-running and recovery

Rerun the same command after interruption. Unchanged completed jobs are reused;
abundance-only changes recalculate expression without remapping proteins.
Use a new `run_name` for fresh work.

ODB retains resumable work in `work/<run_name>/orthogroups/mapping/`;
changed inputs use separate work. Native results/logs also appear in `results/`.
After replacing ODB software in place, use `--forcerun odb_map`.

KofamScan retains work in `work/<run_name>/kegg/`: completed species annotations
are reused, while failed annotations restart with prior attempts retained.
See [reference updates](references.md) when changing snapshots.
