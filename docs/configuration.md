# Configuration

[Documentation](index.md) · [Input formats](inputs.md)

For a new dataset, copy the default file and edit its input paths and OrthoDB node:

```bash
cp config/config.yaml config/mydata.yaml
```

Keep an existing dataset configuration when resuming. `config/mydata.yaml` and
`config/*.local.yaml` are ignored by Git. The complete defaults are in
[config/config.yaml](../config/config.yaml).

## Loading settings and paths

Settings load in this order: `config/config.yaml`, `--configfile` files from left
to right, then command-line `--config` values. Unknown keys are rejected.

```bash
./run_pipeline.sh --configfile config/mydata.yaml config/pilot.yaml \
  --cores 8 --resources mem_gb=64 -- prepare
```

Put options before `--` and targets after it. Paths are relative to the repository
root unless absolute; submit Slurm jobs there or use `sbatch --chdir`.

`run_name` selects directories under `results/`, `work/`, and `logs/`.
Use a new name to retain an earlier analysis. It does not change species selection
or enabled branches. `container_image: auto` selects the matching release image;
see [container setup](containers.md) for overrides.

## Core settings

| Setting | Default | Meaning |
| --- | --- | --- |
| `run_name` | `run001` | Simple directory name: letters, digits, underscores, dots, or hyphens; starts with a letter or digit |
| `container_image` | `auto` | GHCR release: `auto` reads `VERSION`; or a version without `v`, e.g. `"0.1.0"`; see [containers](containers.md) |
| `seed` | `12345` | Shared random seed for supported steps; integer from 1 to 2147483647 |
| `inputs.metadata` | `input/metadata.tsv` | Sample metadata |
| `inputs.busco` | `input/busco/summary.tsv` | BUSCO summary for species selection |
| `inputs.cds_dir` | `input/cds` | Per-species CDS directory |
| `inputs.quant_dir` | `input/quant` | Per-species/per-run abundance directory |
| `inputs.species_trait` | `input/species_trait.tsv` | Species trait table |
| `selection.busco_threshold` | `0.5` | Minimum complete BUSCO fraction |
| `selection.species_list` | `null` | Candidate species file (one ID per line); `null` considers all species. BUSCO filtering follows; see [species selection](inputs.md#species-selection) |
| `selection.missing_taxonomy` | `error` | `error` or `allow` for unresolved taxids |
| `translation.table` | `1` | Genetic code for CDS translation |
| `tpm.multimap` | `error` | `error`, `drop`, or `split` for genes assigned to multiple OGs; see [TPM interpretation](outputs.md#tpm-interpretation) |

## Reproducibility

Top-level `seed` controls VeryFastTree, ASTRAL-IV, representative selection, and
contrast-pair selection, including recomputation by `filter_species`. Each job
initializes its own random generator with this value. The resolved seed is
recorded with the configuration in `run.json` and in the affected steps' reports
or commands. Changing it invalidates affected jobs and their downstream outputs;
unchanged completed jobs remain reusable.

For reproducible results, retain the same inputs, reference snapshots (including
the TimeTree cache), tool versions, actual per-job thread counts, and execution
environment. Seed alone does not ensure identical results across environments.
Logs and reports containing timestamps or temporary paths can differ even when
the scientific results match. Tools without an exposed seed keep their existing
behavior; KofamScan's HMMER uses its fixed default seed.

## OrthoDB settings

Set `odb.node` to a supported OrthoDB v12 mapping level covering your species.
The default `3193` is Embryophyta (land plants). See
[node selection](references.md#choosing-an-orthodb-node) before running a new dataset.

| Setting | Default | Meaning |
| --- | --- | --- |
| `odb.node` | `3193` | NCBI Taxonomy ID of a supported OrthoDB v12 mapping level; Embryophyta by default |
| `odb.existing_results` | `null` | Verified import snapshot of earlier annotations; see [importing ODB results](migration.md#importing-odb-results) |

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

See [targets](running.md#targets) for enabling and requesting analyses, and
[migration](migration.md) for rejected settings in older overrides.
