# Configuration

[Documentation](index.md) · [Input formats](inputs.md)

Edit [config/config.yaml](../config/config.yaml) directly and adjust the settings
as needed. Keep your existing settings when resuming an analysis.

## Loading settings and paths

On Slurm, submit from the repository root:

```bash
sbatch run_pipeline.sh
```

CPU and memory are taken automatically from the Slurm allocation requested by the
script's `#SBATCH` settings. Adjust those settings for your cluster.

`run_name` selects directories under `results/`, `work/`, and `logs/`.
Use a new name to retain an earlier analysis. It does not change species selection
or enabled branches. `container_image: auto` selects the matching release image;
see [container setup](containers.md) for overrides.

## Core settings

| Setting | Default | Meaning |
| --- | --- | --- |
| `run_name` | `run001` | Simple directory name: letters, digits, underscores, dots, or hyphens; starts with a letter or digit |
| `container_image` | `auto` | GHCR release: `auto` reads `VERSION`; or a version without `v`, e.g. `"0.1.0"`; see [containers](containers.md) |
| `inputs.metadata` | `input/metadata.tsv` | Sample metadata |
| `inputs.busco` | `input/busco/summary.tsv` | BUSCO summary for species selection |
| `inputs.cds_dir` | `input/cds` | Per-species CDS directory |
| `inputs.quant_dir` | `input/quant` | Per-species/per-run abundance directory |
| `inputs.species_trait` | `input/species_trait.tsv` | Species trait table |
| `selection.busco_threshold` | `0.5` | Minimum complete BUSCO fraction |

## OrthoDB settings

Set `odb.node` to a supported OrthoDB v12 mapping level covering your species.
The default `3193` is Embryophyta (land plants). See
[node selection](references.md#choosing-an-orthodb-node) before running a new dataset.

| Setting | Default | Meaning |
| --- | --- | --- |
| `odb.node` | `3193` | NCBI Taxonomy ID of a supported OrthoDB v12 mapping level; Embryophyta by default |

CPU and memory defaults live in the rules. Use Snakemake's standard options for
[resource budgets and per-rule overrides](running.md#resource-budgets).
ODB chunks contain up to 100 species. Internal batch size is four times the
allocated ODB thread count; both are defined in [odb.smk](../workflow/rules/odb.smk).
Choose storage and disk capacity for your workload. Work files are retained;
see [re-running and recovery](running.md#re-running-and-recovery).

## Optional analyses and exports

Branch-specific settings are documented with their methods and outputs:

| Configuration section | Guide |
| --- | --- |
| `alignment` | [All-copy OG alignments](alignments.md) |
| `kegg` | [KO annotation and expression](kegg.md) |
| `phylogeny` | [Species sets, markers, sequences, rooting, and resources](phylogeny.md) |
| `phylogeny.dating` | [Manual/TimeTree calibrations and LSD2](dating.md) |
| `contrast` | [Representative analysis and molecular-tree pairs](contrast_pairs.md) |
| `taxonomy_check` | [Taxonomic ranks and MonoPhy review](taxonomy_check.md) |
| `exclude_species` | [Manual exclusion from completed results](species_filter.md) |

The [PhenoRadar input collector](phenoradar_inputs.md) automatically collects
available completed results and needs no configuration section.

See [targets](running.md#targets) for enabling and requesting analyses.
