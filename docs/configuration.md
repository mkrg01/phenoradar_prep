# Configuration

[Documentation](index.md) · [Input formats](inputs.md)

Edit [config/config.yaml](../config/config.yaml) directly. Keep existing settings
when resuming an analysis; see [running the workflow](running.md) for execution.

## Loading settings and paths

Run from the repository root. Dataset paths are fixed under `input/`;
see the [file layout](inputs.md#file-formats).
`config/config.yaml` loads automatically. Supply an override with
`--configfile path/to/override.yaml`; unspecified settings retain their defaults
from the main file. [Optional treePL overrides](dating.md#optional-overrides) use
the workflow's internal defaults when omitted. See the [pilot example](running.md#pilot-run).

`run_name` selects subdirectories under `results/`, `work/`, and `logs/`.
Use a new name to retain an earlier analysis. The container image matches
[VERSION](../VERSION); see [container setup](containers.md).

## Core settings

| Setting | Default | Use |
| --- | --- | --- |
| `run_name` | `run001` | Directory name: letters, digits, underscores, dots, hyphens; starts with a letter/digit |
| `selection.busco_threshold` | `0.5` | Minimum complete BUSCO fraction |
| `selection.species_list` | `false` | Use `input/species_list.txt` when true; [selection rules](inputs.md#species-selection) |
| `translation.table` | `1` | Genetic code for CDS translation |
| `tpm.multimap` | `error` | Ambiguous gene assignments; [TPM policies](outputs.md#tpm-interpretation) |
| `seed` | `12345` | Tree inference, dating, representative selection, and contrast pairs |
| `trait` | `carnivory` | Shared trait column for species selection, pairs, and metadata |

## OrthoDB settings

`odb.node` defaults to `3193` (Embryophyta). Choose a supported OrthoDB v12 level
covering all your species; see [node selection](references.md#choosing-an-orthodb-node).
`odb.existing_results` defaults to `null` (new mapping); set a snapshot directory
to [reuse existing annotations](references.md#reusing-existing-odb-results).
CPU/memory settings use [launcher and rule overrides](running.md#resource-budgets).

## Optional analyses and exports

| Configuration section | Guide |
| --- | --- |
| `alignment` | [All-copy OG alignments](alignments.md) |
| `kegg` | [KO annotation and expression](kegg.md) |
| `phylogeny` | [BUSCO species trees](phylogeny.md) |
| `phylogeny.dating` | [Calibrations and treePL dating](dating.md) |
| `phylogeny.contrast_pairs` | [Trait contrast pairs](contrast_pairs.md) |
| `phylogeny.taxonomy_check` | [MonoPhy review](taxonomy_check.md) |
| `exclude_species` | [Manual species exclusion](species_filter.md) |

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
