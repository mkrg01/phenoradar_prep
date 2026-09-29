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
See [resource budgets](running.md#resource-budgets).

Alternative files are optional: `--config` selects one, and an analysis config
can be a partial override of `config/analysis.yaml`. `WORKFLOW_PYTHON` selects
a host Python with PyYAML. The workflow image follows [VERSION](../VERSION);
GeneGalleon source/image pins belong to build settings.

## Settings by phase

| Build (`build.yaml`) | Analysis (`analysis.yaml`) |
| --- | --- |
| `name` (overridden by `submit --name`), `metadata`, `excluded_accessions`, `store` | Completed `build` or copied `database/` bundle |
| `genegalleon` software and assembly/quant settings | `inputs` for traits, species list, and calibrations |
| `busco.lineage`, `translation.table` | `selection.busco_threshold`, `exclude_species` |
| `odb.node`, `odb.cache_dir`, `odb.chunk_size`, `tpm.multimap` | `trait`, `seed`, optional branches |
| Build `slurm` resources | Analysis `slurm` resources |

Build requires complete products for every included sample. The manual
[accession list](datasets.md#manually-excluding-unusable-accessions) excludes runs
before scheduling; BUSCO acceptance thresholds apply only in analysis.
Downstream inherits lineage, genetic code, ODB node, and TPM ambiguity policy from its database.

Build always disables AMALGKIT rRNA and contamination filtering; these are not
configurable GeneGalleon settings.

Build discovers matching ODB snapshots automatically and maps missing species
in chunks (default 20). `odb.node` defaults to v12 taxid `3193` (Embryophyta).
Retain the configured cache to reuse completed mappings; see [references](references.md).

Analysis uses BUSCO completeness `>= selection.busco_threshold` (default `0.5`).
`selection.species_list: true` enables `inputs.species_list`. Both the list and
`exclude_species` accept biological species IDs (all their samples) or exact
analysis sample IDs, applied before computation. Set `inputs.species_trait: null`
when traits are unused. Build setting `tpm.multimap` defaults to `error`; [OG expression](outputs.md#tpm-interpretation)
also supports `drop` and `split`. Set it in `config/build.yaml` before preparation;
it is not a downstream override. Changing the policy recomputes affected expression
products but reuses assembly, quantification, and mappings. Expression caches are
stored under `<store>/.expression/`, keyed by sample identity, abundance content,
mapping content, policy, and aggregation implementation.

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
