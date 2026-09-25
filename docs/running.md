# Running the workflow

[Documentation](index.md) · [Configuration](configuration.md)

## Installation

Use a published release checkout with the [requirements](../README.md#requirements)
on Linux/Bash, and edit [config/config.yaml](../config/config.yaml).
Place [input files](inputs.md#file-formats) in `input/`.
See [container setup](containers.md) for compute-node tools and images.
First use needs network for the image and missing [references](references.md).
Run `./run_pipeline.sh --prepare-container --cores 1 --resources mem_gb=4`
once before the first container analysis or dry-run; see [container setup](containers.md).
`SNAKEMAKE_BIN` can select another Snakemake executable.

## Slurm

Adjust the `#SBATCH` lines in `run_pipeline.sh` for your cluster (partition,
account, CPUs, memory, time), then submit from the repository root:

```bash
sbatch run_pipeline.sh
```

Monitor `pipeline-<job_id>.out` / `.err`
in the repository root and per-step logs in `logs/<run_name>/`.
Avoid concurrent jobs writing the same results.

## Distributed Slurm execution

The default launcher uses one allocation. To submit workflow rules as independent
Slurm jobs, install `snakemake-executor-plugin-slurm` and use:

```bash
./run_pipeline.sh --slurm --profile profiles/slurm --jobs 10 -- all
```

Run the controller in an appropriate allocated environment. Its resources are
independent of workers; configure worker resources in the profile. Use
[the dataset interface](datasets.md) for GeneGalleon species arrays, frozen input
selection, and staged submission. Do not run copies of the complete workflow as
an array against the same output directory.

## Direct execution

From the repository root, supply a CPU and memory budget:

```bash
./run_pipeline.sh --cores 16 --resources mem_gb=192
```

Put options before `--` and targets after it. Omitting the target runs `all`.
`config/config.yaml` loads automatically; use `--configfile` for
[overrides](configuration.md#loading-settings-and-paths).

## Pilot run

List a few candidate species IDs in `input/species_list.txt`, one per line.
[BUSCO filtering](inputs.md#species-selection) still applies. Then run:

```bash
sbatch --cpus-per-task=8 --mem=80G \
  run_pipeline.sh --configfile config/pilot.yaml \
  --set-threads odb_map=8 --set-resources odb_map:mem_mb=64000
```

[config/pilot.yaml](../config/pilot.yaml) enables that list and sets `run_name: pilot`.
Check selection/mapping QC in `results/pilot/`, runtime, disk use,
and peak memory before a full run. The main configuration leaves the list disabled.

## Targets

Targets schedule missing prerequisites automatically. Phylogeny targets use
`phylogeny.trees` and require the corresponding postprocessing flag for
`contrast_pairs`, `timetree`, `phylogeny_calibrations`, and `taxonomy_check`.
`phylogeny` stops at tree inference; `all` also runs enabled postprocessing.
Alignment and KEGG targets retain their explicit-target opt-in behavior.

| Target | Work requested |
| --- | --- |
| `all` (default) | OG expression/QC and enabled alignment, KEGG, phylogeny, and contrast branches |
| `references`, `kegg_references` | [OrthoDB or KOfam/KEGG snapshots](references.md) |
| `proteins` | CDS translation |
| `mapping` | ODB mapping and merged gene-to-OG index |
| `alignments` | [All-copy OG alignments](alignments.md) |
| `kegg` | [KO annotation and original-TPM sums](kegg.md) |
| `phylogeny_prepare` | BUSCO input audit, outgroup, and marker plan |
| `phylogeny` | Infer the trees selected by `phylogeny.trees` |
| `phylogeny_calibrations` | Trees and [TimeTree calibrations](dating.md#timetree-calibrations), without dating |
| `timetree` | [treePL dating](dating.md), with TimeTree or manual age calibrations |
| `taxonomy_check` | [MonoPhy review](taxonomy_check.md) |
| `contrast_pairs` | [Pairs on the selected trees](contrast_pairs.md); includes their missing inference steps |
| `filter_species` | [Export completed results after exclusions](species_filter.md) |
| `phenoradar_inputs` | [Collect available completed results](phenoradar_inputs.md) |

All phylogeny targets follow `phylogeny.trees`; a target never changes the species
set. Filtering and collection are manual targets that never start producer analyses.

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
The [dating rule](dating.md#resources) requires one thread; its memory can be overridden.
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
To import previously saved annotations, set
[`odb.existing_results`](references.md#reusing-existing-odb-results).

KofamScan retains work in `work/<run_name>/kegg/`: completed species annotations
are reused, while failed annotations restart with prior attempts retained.
See [reference updates](references.md) when changing snapshots.
