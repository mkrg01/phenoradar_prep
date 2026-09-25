# Running the workflow

[Documentation](index.md) · [Configuration](configuration.md)

## Installation and normal execution

Use Linux/Bash with the [requirements](../README.md#requirements) and prepare the
[containers](containers.md). Start with manually curated AMALGKIT metadata and
the [build/analysis commands](datasets.md). Each phase submits its own Slurm
controller; heavy work runs as species arrays or independent rule jobs.

## Direct execution

A prepared analysis can run locally with an explicit budget:

```bash
./run_analysis.sh run --analysis analyses/analysis001 --local --cores 16 --mem-mb 192000
```

`run_pipeline.sh` remains a low-level launcher. Supply the **resolved** phase
configuration, not `config/build.yaml` or `config/analysis.yaml`:

```bash
./run_pipeline.sh --cores 16 --resources mem_gb=192 \
  --configfile analyses/analysis001/pipeline.yaml -- all
```

Put options before `--` and targets after it. A resolved analysis configuration
can only consume completed build products. A resolved build configuration is
used for the `mapping` target; afterwards publish its verified completion with
`run_build.sh complete --build builds/<id>`. Normal build submission does this
automatically. For an existing schema-1 dataset use its original checkout to
resume, or [migrate](datasets.md#migrating-old-configurations).

## Pilot run

For upstream inspection, pass `--species-list pilot.txt` to build `submit`, then
resume without that list to finish all metadata species. For a small downstream
analysis, [config/pilot.yaml](../config/pilot.yaml) is an analysis override:

```bash
./run_analysis.sh prepare --config config/pilot.yaml \
  --build builds/build001 --name pilot
./run_analysis.sh submit --analysis analyses/pilot
```

Prepare `input/species_list.txt` with candidate species IDs first; the BUSCO
threshold still applies. Inspect selection and mapping QC before larger analyses.

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
Rule defaults can be overridden by the relevant phase's `slurm.rules`; `runtime` is in minutes.

| Rule | Job unit | Threads | Memory (GB) |
| --- | --- | --- | --- |
| `odb_map` | Mapping chunk (default 20 species) | 16 | 192 |
| `align_orthogroup` | OG | 4 | 8 |
| `annotate_kofam` | Species | 4 | 8 |
| `check_taxonomy` | Species set | 1 | 8 |

See [phylogeny resources](phylogeny.md#resources) for tree-inference defaults.
The [dating rule](dating.md#resources) requires one thread; its memory can be overridden.
Two default ODB jobs need 32 CPUs and 384 GB. To reduce per-job requests:

```bash
./run_pipeline.sh --cores 16 --resources mem_gb=192 \
  --configfile builds/build001/pipeline.yaml \
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
Use a new build or analysis ID for changed scientific conditions.

ODB retains resumable work in `work/<run_name>/orthogroups/mapping/`;
changed inputs use separate work. Native results/logs also appear in `results/`.
Use a separate build/store/cache for deliberate changes to upstream software; do not edit frozen runs.
To import previously saved annotations, set
[`odb.existing_results`](references.md#reusing-existing-odb-results).

KofamScan retains work in `work/<run_name>/kegg/`: completed species annotations
are reused, while failed annotations restart with prior attempts retained.
See [reference updates](references.md) when changing snapshots.
