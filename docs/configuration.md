# Configuration

[Documentation](index.md) · [Input formats](inputs.md)

Edit [build.yaml](../config/build.yaml) for reusable sample products and
[analysis.yaml](../config/analysis.yaml) for downstream analyses. Both commands
load these files by default. See the [build/analysis guide](datasets.md) for commands.

## Loading settings and paths

Configuration paths are relative to the repository root; local FASTQ paths in
metadata are relative to the metadata file's directory. The first `submit` freezes
settings and inputs under `results/<build>/` or `results/<build>/downstream/<analysis>/`.
Submitting the same name reuses those conditions, including after `--dry-run`.
Later source edits affect new runs only; scientific changes require a new name. Do not edit generated
`pipeline.yaml` files. `workflow/pipeline_defaults.yaml` is an internal schema.

For retries with different CPU, memory, time, or concurrency, edit the phase
config and pass `submit --resources config/build.yaml` or
`submit --resources config/analysis.yaml`. Only its `slurm` section is applied.
Use `slurm.total_limits` for the job count and total CPU/memory budgets of one
build or analysis, and `slurm.per_job_resources.<job>` for individual job requests.
Total CPU/memory budgets include the controller; `null` adds no workflow-wide
limit. The job-count budget also defaults to `null`. Unspecified worker settings use internal defaults, so a separate
`default_resources` section is unnecessary.
Per-job memory uses `mem_gb` (positive integer GB), for example `mem_gb: 128`
in place of `mem_mb: 128000`. Job durations use a quoted Slurm `time`, for example
`time: "1-00:00:00"` in place of `runtime: 1440`.
See [resource budgets](running.md#resource-budgets).

Alternative files are optional: `--config` selects one, and an analysis config
can be a partial override of `config/analysis.yaml`. `WORKFLOW_PYTHON` selects
a host Python with PyYAML. The workflow image follows [VERSION](../VERSION);
GeneGalleon source/image pins belong to build settings.
`genegalleon.image_uri` must include an OCI SHA256 digest in the form
`docker://<registry>/<image>@sha256:<64 lowercase hexadecimal characters>`;
omitting it, using a tag alone, or using a direct HTTPS SIF URL is an error.
The separate `image_sha256` setting is no longer accepted. See
[GeneGalleon software retrieval and caching](datasets.md#pinned-genegalleon-source-and-sif).

## Settings by phase

| Build (`build.yaml`) | Analysis (`analysis.yaml`) |
| --- | --- |
| `name` (overridden by `submit --name`), `metadata`, `excluded_accessions`, `reuse_from` | Completed `build` or copied `database/` bundle |
| `genegalleon` software and assembly/quant settings | `inputs` for traits, species list, and calibrations |
| `busco.lineage`, `translation.table` | `selection.busco_threshold`, `exclude_species` |
| `odb.ncbi_tax_id`, `odb.chunk_size` | `trait`, `seed`, optional branches |
| Build `slurm` resources | Analysis `slurm` resources |

Build requires complete products for every included sample. The manual
[accession list](datasets.md#manually-excluding-unusable-accessions) excludes runs
before scheduling; BUSCO acceptance thresholds apply only in analysis.
Downstream inherits lineage, genetic code, and the ODB mapping NCBI Taxonomy ID from its database.
Use `odb.ncbi_tax_id` in configs; this replaces the former `odb.node` key.
Existing database/reference snapshots retain their `node` metadata field.

Build always disables AMALGKIT rRNA and contamination filtering; these are not
configurable GeneGalleon settings.

Build reuses matching sample products from `reuse_from`, which accepts `null`,
a completed database path, or a list of paths. Paths are relative to the project
root; `results/<build>/database/` and its parent build directory are accepted.
See [database reuse](datasets.md#reusing-completed-databases) for examples and conflicts.
Missing samples are computed; ODB mapping runs in chunks (default 50), using
v12 taxid `3193` (Embryophyta) by default. The build manages its working caches
under `results/<name>/work/cache/`. Separate `store` and `odb.cache_dir` settings
are no longer accepted in `build.yaml`.

Analysis uses BUSCO completeness `>= selection.busco_threshold` (default `0.5`).
`selection.species_list: true` enables `inputs.species_list`. Both the list and
`exclude_species` accept biological species IDs (all their samples) or exact
analysis sample IDs, applied before computation. Set `inputs.species_trait: null`
when traits are unused. [OG expression](outputs.md#tpm-interpretation) always stops
with an error if a quantified gene maps to multiple OGs; this behavior is not
configurable. Expression caches are stored under `results/<build>/work/cache/expression/`, keyed by
sample identity, abundance content, mapping content, and aggregation implementation.

## Optional analyses and exports

Set options **before the first analysis submission**, including `--dry-run`. Target examples in the guides
assume an already prepared `results/angiosperm_leaf_20260925/downstream/c4_photosynthesis_20260929`.

| Configuration section | Guide |
| --- | --- |
| `alignment` | [All-copy OG alignments](alignments.md) |
| `kegg` | [KO annotation and expression](kegg.md) |
| `phylogeny` | [BUSCO species trees](phylogeny.md) |
| `phylogeny.dating` | [Calibrations and dating](dating.md) |
| `phylogeny.contrast_pairs` | [Trait contrast pairs](contrast_pairs.md) |
| `phylogeny.taxonomy_check` | [Taxonomic review](taxonomy_check.md) |

[PhenoRadar collection](phenoradar_inputs.md) needs no settings.
[Post hoc species filtering](species_filter.md) is a separate export utility.

## Choosing species trees

`phylogeny.trees` selects `all`, `phenotyped`, `representatives`, or several trees
for independent inference; `[]` disables inference. See
[species sets](phylogeny.md#species-sets) for selection details.

Enabled `contrast_pairs`, `dating`, and `taxonomy_check` postprocessing applies to
those trees and requires a nonempty list. Dating and taxonomy checks support only
`all` and `phenotyped`. Explicit postprocessing targets require their enabled flag;
`phylogeny` stops at inference, while `all` includes enabled postprocessing.
