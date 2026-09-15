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

## OrthoDB settings

`odb.node` defaults to `3193` (Embryophyta). Choose a supported OrthoDB v12 level
covering all your species; see [node selection](references.md#choosing-an-orthodb-node).
CPU/memory settings use [launcher and rule overrides](running.md#resource-budgets).

## Optional analyses and exports

| Configuration section | Guide |
| --- | --- |
| `alignment` | [All-copy OG alignments](alignments.md) |
| `kegg` | [KO annotation and expression](kegg.md) |
| `phylogeny` | [BUSCO species trees](phylogeny.md) |
| `phylogeny.dating` | [Calibrations and treePL dating](dating.md) |
| `contrast` | [Trait contrast pairs](contrast_pairs.md) |
| `taxonomy_check` | [MonoPhy review](taxonomy_check.md) |
| `exclude_species` | [Manual species exclusion](species_filter.md) |

See [targets](running.md#targets) for requesting analyses.
[PhenoRadar input collection](phenoradar_inputs.md) needs no configuration section.
