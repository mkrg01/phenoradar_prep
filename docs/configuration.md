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

The workflow loads `config/config.yaml`, then overlays files supplied with
`--configfile` in order. Command-line `--config` values take precedence.
Unknown keys, including retired settings, are rejected before workflow execution
with their full configuration path. Configuration sections must be mappings.

```bash
./run_pipeline.sh --software-deployment-method conda \
  --configfile config/mydata.yaml config/pilot.yaml \
  --cores 8 --resources mem_gb=64 -- prepare
```

Put all options before `--` and targets after it. `--configfile` accepts multiple
filenames, so a target placed next to a filename can be read as another file.

`analysis` selects a directory name under `results/`, `work/`, and `logs/`.
Use a new name to retain results for comparison. The default is `full`;
[config/pilot.yaml](../config/pilot.yaml) uses `pilot` and a user-supplied species
list. The name itself does not change species selection.

Paths are relative to the repository root unless absolute. Direct execution of
`run_pipeline.sh` changes to that directory. Submit Slurm jobs from the repository
root, or use `sbatch --chdir=/path/to/phenoradar_prep`. Generated storage locations
are fixed; see the [directory layout](outputs.md#directory-layout).

## Core settings

| Setting | Default | Meaning |
| --- | --- | --- |
| `analysis` | `full` | Simple directory name: letters, digits, underscores, dots, or hyphens; starts with a letter or digit |
| `container_image` | `null` | Image URI or absolute SIF path for `--sdm conda apptainer`; see [containers](containers.md) |
| `inputs.metadata` | `input/metadata.tsv` | Sample metadata |
| `inputs.busco` | `input/busco/summary.tsv` | BUSCO summary for species selection |
| `inputs.cds_dir` | `input/cds` | Per-species CDS directory |
| `inputs.quant_dir` | `input/quant` | Per-species/per-run abundance directory |
| `inputs.species_trait` | `input/species_trait.tsv` | Species trait table |
| `selection.busco_threshold` | `0.5` | Minimum complete BUSCO fraction |
| `selection.species_list` | `null` | Optional list restricting the eligible species |
| `selection.missing_taxonomy` | `error` | `error` or `allow` for unresolved taxids |
| `taxonomy.source` | `null` | Local ETE4-compatible SQLite source; otherwise download NCBI taxonomy when the snapshot is missing |
| `translation.table` | `1` | Genetic code for CDS translation |
| `tpm.multimap` | `error` | `error`, `drop`, or `split` for genes assigned to multiple OGs; see [TPM interpretation](outputs.md#tpm-interpretation) |

## OrthoDB settings

The supported OrthoDB release is v12. **Review `odb.node` for your dataset**;
the default `3193` is dataset-specific and selects `resources/orthodb/v12_3193/`.

| Setting | Default | Meaning |
| --- | --- | --- |
| `odb.node` | `3193` | Positive integer taxonomic node |
| `odb.chunk_size` | `50` | Maximum species per mapping chunk |
| `odb.threads` | `16` | Workers per chunk |
| `odb.batch_size` | `64` | Batch size; must be at least `odb.threads` |
| `odb.mem_gb` | `192` | Scheduling memory per chunk, in decimal GB |
| `odb.min_free_gb` | `750` | Required free disk space before mapping, in GiB |
| `odb.reference_min_free_gb` | `200` | Required free disk space before reference preparation, in GiB |
| `odb.allow_nonlocal` | `false` | Permit mapping work on filesystems outside the supported local types |
| `odb.keep_work` | `false` | Retain successful mapping work directories |
| `odb.existing_results` | `null` | Verified import snapshot of earlier annotations; see [importing ODB results](migration.md#importing-odb-results) |

Memory settings use positive whole decimal GB (1 GB = 1000 MB). See
[resource budgets](running.md#resource-budgets) for allocation sizing and concurrency.

## Optional analyses and exports

Branch-specific settings are documented with their methods and outputs:

| Configuration section | Guide |
| --- | --- |
| `alignment` | [All-copy OG alignments](alignments.md) |
| `kegg` | [KO annotation and expression](kegg.md) |
| `phylogeny` | [Species sets, markers, sequences, rooting, and resources](phylogeny.md) |
| `phylogeny.dating` | [Manual/TimeTree calibrations and LSD2](dating.md) |
| `contrast` | [Representative analysis and molecular-tree pairs](contrast_pairs.md) |
| `taxonomy_audit` | [Taxonomic ranks and MonoPhy review](taxonomy_audit.md) |
| `exclude_species` | [Manual exclusion from completed results](species_filter.md) |
| `phenoradar` | [Selection of completed downstream inputs](phenoradar_inputs.md) |

See [targets](running.md#targets) for enabling and requesting analyses, and
[migration](migration.md) for rejected settings in older overrides.
