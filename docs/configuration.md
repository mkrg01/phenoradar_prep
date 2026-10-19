# Configuration

[Documentation](index.md) · [Input formats](inputs.md)

Edit [config/build.yaml](../config/build.yaml) for sample processing and
[config/analysis.yaml](../config/analysis.yaml) for downstream analyses.
`--config PATH` selects another file. Analysis overrides can contain only the
settings you want to change.

## Loading settings and paths

Config paths are relative to the repository root. Local FASTQ paths are relative
to the metadata file's directory. Keep imported databases inside the repository
so the container can access them.

The first submission saves inputs and settings, including with `--dry-run`.
Use a new run for changed scientific conditions; use
[resource overrides](running.md#resource-budgets) for job resources.

## Build settings

| Setting | Purpose |
| --- | --- |
| `name` | Dataset prefix; exact directory name with `name_mode: fixed`. `--name` overrides it |
| `name_mode` | `timestamp` (supplied config) creates dated builds; `fixed` resumes the named build |
| `metadata`, `excluded_accessions` | Sample metadata and optional run exclusions |
| `metadata_provenance` | Optional metadata history JSON, matched by checksum and copied into the database |
| `reuse_from` | Previous build/database path, list of paths, `auto`, or `null`; see [reuse](datasets.md#reusing-completed-databases) |
| `rebuild_samples` | Sample IDs or run accessions to restart from assembly in a new build; default `[]` |
| `genegalleon` | [Pinned software](datasets.md#pinned-genegalleon-source-and-sif) and assembly/quantification settings |
| `busco.lineage` | BUSCO dataset; default `embryophyta_odb12` |
| `busco.cache_dir` | Shared reference directory; default `resources/busco`; `null` disables sharing |
| `translation.table` | NCBI genetic code; default `1` |
| `slurm` | [Job resources and total limits](running.md#resource-budgets) |

Resume saved builds with `submit --build PATH`. See
[naming and previews](datasets.md#saved-settings-and-previews).

### OrthoDB

All samples use the selected release and clade. In a **new build**, samples already
processed with that reference are reused; the others are mapped automatically.

| Setting | Default | Purpose |
| --- | --- | --- |
| `odb.version` | `v12` | Fixed release supported by the installed ODB-mapper |
| `odb.ncbi_tax_id` | `3193` | [Mapping clade](references.md#choosing-an-orthodb-mapping-clade); `3193` is Embryophyta |
| `odb.chunk_size` | `50` | Maximum samples per mapping chunk |

See [OrthoDB update examples](datasets.md#update-orthodb-for-all-samples).

### Shared BUSCO references

New builds download each lineage once to `resources/busco/<lineage>/` and
share it read-only across samples. This requires Apptainer; for Docker or
SingularityCE, set `busco.cache_dir: null` to use per-sample downloads.

References are reused without automatic updates and checked for changes.
Keep the cache while processing or resuming builds. To download a different
snapshot, choose a new `busco.cache_dir` in a new build; use `rebuild_samples`
to recompute completed samples.

Existing per-sample copies are not moved or deleted.

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

BUSCO lineage, genetic code, and OrthoDB reference are inherited from the build.
Set unused input paths to `null` and disable branches that require those inputs.

## Intermediate storage

Both build and analysis configs use `storage.keep_intermediates: false` by
default, enabling [cleanup](outputs.md#storage-cleanup). GeneGalleon computation scratch
is removed after successful and failed sample jobs; native restart checkpoints
are preserved.
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
