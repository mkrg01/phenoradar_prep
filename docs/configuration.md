# Configuration

[Documentation](index.md) · [Input formats](inputs.md)

Edit [build.yaml](../config/build.yaml) for reusable species products and
[analysis.yaml](../config/analysis.yaml) for downstream analyses. Both commands
load these files by default. See the [build/analysis guide](datasets.md) for commands.

## Loading settings and paths

Paths are relative to the repository root. `prepare` freezes settings and inputs
under `builds/<build>/` or `analyses/<analysis>/`. Later source edits affect new
preparations only; scientific changes require a new name. Do not edit generated
`pipeline.yaml` files. `workflow/pipeline_defaults.yaml` is an internal schema.

For retries with different CPU, memory, time, or concurrency, edit the phase
config and pass `submit --resources config/build.yaml` or
`submit --resources config/analysis.yaml`. Only its `slurm` section is applied.
See [resource budgets](running.md#resource-budgets).

Alternative files are optional: `--config` selects one, and an analysis config
can be a partial override of `config/analysis.yaml`. `WORKFLOW_PYTHON` selects
a host Python with PyYAML. The workflow image follows [VERSION](../VERSION);
GeneGalleon source/image pins belong to build settings.

## Settings by phase

| Build (`build.yaml`) | Analysis (`analysis.yaml`) |
| --- | --- |
| `name` (overridden by `prepare --name`), `metadata`, `excluded_accessions`, `store` | Completed `build` or copied `products/` bundle |
| `genegalleon` software and assembly/quant settings | `inputs` for traits, species list, and calibrations |
| `busco.lineage`, `translation.table` | `selection.busco_threshold`, `exclude_species` |
| `odb.node`, `odb.cache_dir`, `odb.chunk_size` | `trait`, `seed`, `tpm`, optional branches |
| Build `slurm` resources | Analysis `slurm` resources |

Build requires complete products for every included species. The manual
[accession list](datasets.md#manually-excluding-unusable-accessions) excludes runs
before scheduling; BUSCO acceptance thresholds apply only in analysis.
Analysis inherits lineage, genetic code, and ODB node from its build.

Build always disables AMALGKIT rRNA and contamination filtering; these are not
configurable GeneGalleon settings.

Build discovers matching ODB snapshots automatically and maps missing species
in chunks (default 20). `odb.node` defaults to v12 taxid `3193` (Embryophyta).
Import old mappings with `register --odb-results <snapshot> --odb-only`;
see [references](references.md).

Analysis uses BUSCO completeness `>= selection.busco_threshold` (default `0.5`).
`selection.species_list: true` enables `inputs.species_list`; `exclude_species`
removes exact IDs before computation. Set `inputs.species_trait: null` when traits
are unused. `tpm.multimap` defaults to `error`; [OG expression](outputs.md#tpm-interpretation)
also supports `drop` and `split`.

## Optional analyses and exports

Set options **before preparing the analysis**. Target examples in the guides
assume an already prepared `analyses/analysis001`.

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
