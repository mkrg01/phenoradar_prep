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
| `container_image` | `auto` | Image matching `VERSION`; override with a release URI or absolute SIF path; see [containers](containers.md) |
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

Set `odb.node` to a supported OrthoDB v12 mapping level covering your species.
The default `3193` is Embryophyta (land plants). See
[node selection](references.md#choosing-an-orthodb-node) before running a new dataset.

| Setting | Default | Meaning |
| --- | --- | --- |
| `odb.node` | `3193` | NCBI Taxonomy ID of a supported OrthoDB v12 mapping level; Embryophyta by default |
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
