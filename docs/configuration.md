# Configuration

[Documentation](index.md) · [Input formats](inputs.md)

Edit [build.yaml](../config/build.yaml) for reusable species products and
[analysis.yaml](../config/analysis.yaml) for analyses of a completed build.
Keep local copies as `config/build.local.yaml` and `config/analysis.local.yaml`.
See the [two-phase guide](datasets.md) for preparation, submission, and migration.

## Loading settings and paths

Paths are relative to the repository root. Build settings are frozen by
`run_build.sh prepare --name <build>`. Analysis settings and auxiliary inputs
are frozen by `run_analysis.sh prepare --name <analysis>`; an analysis config can
be a partial override of `config/analysis.yaml`. Execution uses generated
`builds/<id>/pipeline.yaml` or `analyses/<id>/pipeline.yaml` files.
`workflow/pipeline_defaults.yaml` is an internal compatibility schema for
Snakemake, not a third user configuration.

Scientific changes require a new build or analysis ID. Retry CPU, memory, time,
and concurrency changes use `submit --resources <yaml>` and are recorded per
submission. The workflow container follows [VERSION](../VERSION);
GeneGalleon source/SIF pins belong to build settings.

## Settings by phase

| Build settings | Analysis settings |
| --- | --- |
| `metadata`, `excluded_accessions`, `store` | Completed `build` |
| `genegalleon` source, SIF and assembly/quant settings | `inputs` for traits, species list and calibrations |
| `busco.lineage`, `translation.table` | `selection.busco_threshold`, `exclude_species` |
| `odb.node`, `existing_results`, `cache_dir`, `chunk_size` | `trait`, `seed`, `tpm`, optional analysis branches |
| Build `slurm` resources | Analysis `slurm` resources |

Build always computes full BUSCO tables for every included metadata species.
`excluded_accessions` applies manual run-level exclusions before scheduling;
see [failure recovery and exclusions](datasets.md#failure-and-recovery-behavior). Acceptance
thresholds apply only in analysis. Genetic code and lineage are inherited from
the completed build and cannot be overridden by analysis.

`odb.node` defaults to OrthoDB v12 taxid `3193` (Embryophyta). Build always uses
incremental mapping: matching imported/native snapshots are reused and missing
species are mapped in chunks (default 20). See [reference configuration](references.md).

`selection.busco_threshold` defaults to `0.5`. `selection.species_list: true`
enables `inputs.species_list`. `exclude_species` removes listed IDs before any
analysis. Set `inputs.species_trait: null` when no trait table is needed;
trait-based tree settings require an actual table. `tpm.multimap` remains
`error` by default; `drop` and `split` are also supported.

## Optional analyses and exports

| Configuration section | Guide |
| --- | --- |
| `alignment` | [All-copy OG alignments](alignments.md) |
| `kegg` | [KO annotation and expression](kegg.md) |
| `phylogeny` | [BUSCO species trees](phylogeny.md) |
| `phylogeny.dating` | [Calibrations and treePL dating](dating.md) |
| `phylogeny.contrast_pairs` | [Trait contrast pairs](contrast_pairs.md) |
| `phylogeny.taxonomy_check` | [MonoPhy review](taxonomy_check.md) |
| `exclude_species` | [Species selection](datasets.md#run-an-analysis) |

See [targets](running.md#targets) for requesting analyses.
[PhenoRadar input collection](phenoradar_inputs.md) needs no configuration section.

## Choosing species trees

`phylogeny.trees` is the single selection of trees to infer. An empty list disables
inference; each selected tree is included in `all` and in the `phylogeny` target.
The common `trait` is read from `input/species_trait.tsv`.

| `phylogeny.trees` | Inference species |
| --- | --- |
| `[]` | No species trees |
| `[all]` | Every species passing input selection and BUSCO filtering |
| `[phenotyped]` | Selected species with a known `trait`, including state zero |
| `[representatives]` | Representatives of homogeneous trait clades on the NCBI guide |
| `[all, representatives]` | Two independent inference runs |

For a compressed tree with contrast pairs:

```yaml
trait: carnivory
phylogeny:
  trees: [representatives]
  contrast_pairs:
    enabled: true
  dating:
    enabled: false
  taxonomy_check:
    enabled: false
```

Representatives are selected from known-trait species; both trait states are
compressed. This selection does not preserve every positive-trait species.
The remaining BUSCO, alignment, rooting, and inference settings apply to every
selected tree. Changing which trees are requested does not change inference
parameters for an existing tree.

`phylogeny.contrast_pairs.enabled` selects pairs on those trees; it never chooses
a different inference species set. `phylogeny.dating.enabled` and
`phylogeny.taxonomy_check.enabled` likewise request postprocessing of each tree.
Dating and taxonomy checks currently support only `all` and `phenotyped`;
combining either with `representatives` is rejected.

Enabled postprocessing requires a nonempty `trees` list. Explicit phylogeny targets
respect these settings: a disabled step produces an error rather than becoming
enabled by its target name. `phylogeny` runs through tree inference only; `all`
also runs enabled postprocessing. Startup logs show the selected trees and steps,
and resolved marker plans report the number of species before inference.
