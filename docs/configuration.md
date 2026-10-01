# Configuration

[Documentation](index.md) · [Input formats](inputs.md)

Edit [config/build.yaml](../config/build.yaml) for reusable sample products and
[config/analysis.yaml](../config/analysis.yaml) for downstream analyses. The
launchers load these files by default; `--config` selects another file. Analysis
configs may partially override the default analysis config.

## Loading settings and paths

Config paths are relative to the repository root. Local FASTQ paths are relative
to the metadata file's directory. Keep imported databases inside the repository
so the container can access them.

The first submission saves inputs and settings, including with `--dry-run`.
Use a new name for changed scientific conditions. See
[saved settings](datasets.md#saved-settings-and-previews) and
[resource overrides](running.md#resource-budgets).

## Build settings

| Setting | Purpose |
| --- | --- |
| `name` | Build directory under `results/`; `submit --name` overrides it |
| `metadata`, `excluded_accessions` | Sample metadata and optional run exclusions |
| `reuse_from` | Completed database path, list of paths, or `null` for a fresh build |
| `genegalleon` | [Pinned software](datasets.md#pinned-genegalleon-source-and-sif) and assembly/quantification settings |
| `busco.lineage` | BUSCO dataset; default `embryophyta_odb12` |
| `translation.table` | NCBI genetic code; default `1` |
| `odb.ncbi_tax_id` | [OrthoDB mapping clade](references.md#choosing-an-orthodb-mapping-clade); default `3193` (Embryophyta) |
| `odb.chunk_size` | Maximum samples per mapping chunk; default `50` |
| `slurm` | [Job resources and total limits](running.md#resource-budgets) |

## Analysis settings

| Setting | Purpose |
| --- | --- |
| `build` | Completed build or copied database; `--build` overrides it |
| `inputs.species_trait`, `trait` | Trait table and column; set the path to `null` when unused |
| `inputs.species_list` | Restrict to listed species/sample IDs; `null` disables filtering |
| `exclude_species` | Omit biological species IDs (all their samples) or exact sample IDs |
| `selection.busco_threshold` | Minimum complete BUSCO fraction per sample; default `0.5` |
| `inputs.calibrations` | Age bounds when dating uses `calibration_source: file` |
| `slurm` | Job resources and total limits for this analysis |

Lineage, genetic code, and ODB clade are inherited from the build.
Configured input paths must exist; set unused optional paths to `null`.
Trait-dependent branches also need to be disabled when traits are unused.

## Intermediate storage

Both build and analysis configs use `storage.keep_intermediates: false` by
default, enabling [cleanup after successful jobs](outputs.md#storage-cleanup).
Set it to `true` before preparing a run to retain intermediates for debugging.

## Optional analyses and exports

Choose branches before the first submission. The supplied analysis config enables
representative trees and contrast pairs; review these choices for your dataset.

| Setting | Guide |
| --- | --- |
| `alignment.enabled` | [All-copy OG alignments](alignments.md) |
| `kegg.enabled`, `kegg.ambiguity` | [KO annotation and expression](kegg.md) |
| `phylogeny.trees`, marker and outgroup settings | [BUSCO trees](phylogeny.md) |
| `phylogeny.dating` | [Calibrations and dating](dating.md) |
| `phylogeny.contrast_pairs` | [Trait pairs](contrast_pairs.md) |
| `phylogeny.taxonomy_check` | [Taxonomic review](taxonomy_check.md) |

`all` runs enabled branches and collects [PhenoRadar inputs](phenoradar_inputs.md).
See [targets](running.md#targets) to run branches separately.

## Choosing species trees

`phylogeny.trees` accepts `all`, `phenotyped`, `representatives`, or several of
these for independent inference. Use `[]` to disable inference and disable any
postprocessing flags as well. Dating and taxonomy checks support only `all` and
`phenotyped`; contrast pairs support all three.

See [species sets](phylogeny.md#species-sets) for selection criteria. Representative
selection, pairs, and tree inference use a fixed internal seed (`12345`).
