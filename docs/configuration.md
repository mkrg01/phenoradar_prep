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
| `taxonomy.cache_dir` | Shared ETE/NCBI database; default `resources/taxonomy`; `null` disables sharing |
| `translation.table` | NCBI genetic code; default `1` |
| `execution` | [Job resources and total limits](running.md#resource-budgets) |

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
share it read-only across samples. Sharing requires Apptainer; for Docker or
SingularityCE, set both `busco.cache_dir` and `taxonomy.cache_dir` to `null`.

References are reused without automatic updates and checked for changes.
Keep the cache while processing or resuming builds. To download a different
snapshot, choose a new `busco.cache_dir` in a new build; use `rebuild_samples`
to recompute completed samples.

### Shared taxonomy database

New builds prepare one ETE/NCBI database under `resources/taxonomy/ete/` and
share it read-only across samples. This requires Apptainer; set
`taxonomy.cache_dir: null` to use per-sample downloads. The snapshot is verified
and reused without automatic updates. Keep it available while processing or
resuming samples; choose a new cache directory for a different snapshot.

## Analysis settings

| Setting | Purpose |
| --- | --- |
| `build` | Completed build or copied database; `--build` overrides it |
| `inputs.species_trait`, `trait` | Trait table and column; set the path to `null` when unused |
| `trait_encoding` | Explicit text-state mapping to integers `0` and `1`; `null` for existing `0`/`1` input |
| `inputs.species_list` | Restrict to listed species/sample IDs; `null` disables filtering |
| `exclude_species` | Omit biological species IDs (all their samples) or exact sample IDs |
| `selection.busco_threshold` | Minimum complete BUSCO fraction per sample; default `0.5` |
| `inputs.calibrations` | Age bounds when dating uses `calibration_source: file` |
| `execution` | Job resources and total limits for this analysis |

BUSCO lineage, genetic code, and OrthoDB reference are inherited from the build.
Set unused input paths to `null` and disable branches that require those inputs.

## Intermediate storage

Both build and analysis configs use `storage.keep_intermediates: false` by
default, enabling [cleanup](outputs.md#storage-cleanup). GeneGalleon computation scratch
is removed after successful and failed sample jobs; unfinished samples retain
their restart data. Completed samples also retire [unneeded work](outputs.md#storage-cleanup).
Set it to `true` before preparing a run to retain intermediates for debugging.

The following policies apply to both schedulers:

| Setting | Default | Purpose |
| --- | --- | --- |
| `storage.scratch_root` | `auto` | Choose `/data1`, node `$TMPDIR`, or `/tmp` on a different filesystem. `env`, `workspace`, or an absolute directory can override it |
| `storage.allow_shared_scratch` | `true` | Permit the automatic project-scratch fallback |
| `storage.compact_logs` | `true` | Archive ODB RunLogs per chunk while retaining their checksums |
| `storage.max_bytes`, `storage.max_inodes` | `null` | Optional project-wide allocated-byte and inode ceilings |
| `storage.interval_seconds` | `60` | Storage inventory interval |
| `storage.scratch_min_free_bytes` | `0` | Required free bytes when selecting scratch and while running |
| `storage.scratch_min_free_inodes` | `0` | Required free inodes, checked before starting and at least once per second between inventories |
| `storage.scratch_max_bytes` | `null` | Per-job scratch ceiling |
| `storage.scratch_max_inodes` | `null` | Per-job scratch inode ceiling |
| `storage.scratch_reserve_bytes`, `storage.scratch_reserve_inodes` | `0` | Reserve headroom against other cooperating jobs on the same node/filesystem before starting |
| `storage.start_fraction`, `storage.stop_fraction` | `0.85`, `0.95` | Refuse new work and stop active work before the configured project ceilings |
| `storage.series_max_bytes` | `1048576` | Bound each live measurement series; keep one compressed previous window and compress on completion |
| `storage.task_limits` | `{}` | Override scratch limits/reservations by task or rule name, e.g. `sample` and `odb_map`, on both schedulers |
| `storage.projection_factor` | `1.3` | Safety factor for pilot extrapolation |

Scratch is selected on each execution node and cleaned after the attempt.
Restart products remain on shared storage. Keeping intermediates forces project
scratch and retains it. Interrupted jobs killed without a cleanup opportunity can
leave scratch directories; their `owner.json` identifies the run.

Image pulls also use the supplied task scratch or temporary directory for unpacking
and disposable OCI layers. Only the completed SIF is published to the image cache.
Explicit `APPTAINER_CACHEDIR`/`SINGULARITY_CACHEDIR` locations remain caller-owned.

The supplied configs, including angiosperm, measure storage without reserving
capacity or enforcing byte/inode ceilings. Normal cleanup and execution locks
remain enabled. To use shared computation storage directly, set
`storage.scratch_root: workspace` or an existing absolute directory; `auto`
prefers node-local storage and permits a project-scratch fallback.

Storage limits are optional. If configured, project ceilings cover the repository,
including other builds and shared caches. New jobs refuse to start at 85% of a
ceiling; active jobs stop at an observed 95%. Reservations coordinate this user's
jobs on each node/filesystem, count in full against current free space, and do
not reserve space against unrelated applications or users.

To retry with a different scratch location or budget, set
`PHENORADAR_STORAGE_POLICY=/absolute/path/storage-retry.yaml` when submitting.
The file must contain only `storage:` settings; `task_limits` entries merge by name.
This changes the operational policy without editing frozen build inputs. Each
usage receipt records the policy checksum and effective limits. To remove the
capacity limits saved by an older build while retaining measurement and cleanup:

```bash
PHENORADAR_STORAGE_POLICY=config/storage-unrestricted.yaml \
  ./run_build.sh submit --build results/{saved_build}
```

Saved CPU, memory, time and concurrency requests are unchanged by that policy;
use [resource overrides](running.md#resource-budgets) to change them.
See [pilot measurements](outputs.md#storage-pilot).

## Optional analyses and exports

Choose branches before the first submission. The supplied analysis config enables
representative trees and contrast pairs; review these choices for your dataset.

| Setting | Guide |
| --- | --- |
| `alignment.enabled`, `alignment.min_sample_fraction` | [All-copy OG alignments](alignments.md); default sample occupancy 0.5 |
| `kegg.enabled`, `kegg.ambiguity` | [KO annotation and expression](kegg.md) |
| `phylogeny.trees`, marker and outgroup settings | [BUSCO trees](phylogeny.md) |
| `phylogeny.dating` | [Calibrations and dating](dating.md) |
| `phylogeny.contrast_pairs` | [Trait pairs](contrast_pairs.md) |
| `phylogeny.taxonomy_check` | [Taxonomic review](taxonomy_check.md) |

`all` runs enabled branches and validates their outputs.
[Use native outputs directly in PhenoRadar](phenoradar_inputs.md).
See [targets](running.md#targets) to run branches separately.

## Choosing species trees

`phylogeny.trees` accepts `all`, `phenotyped`, `representatives`, or several of
these for independent inference. Use `[]` to disable inference and disable any
postprocessing flags as well. Dating and taxonomy checks support only `all` and
`phenotyped`; contrast pairs support all three.

See [species sets](phylogeny.md#species-sets) for selection criteria. Representative
selection, pairs, and tree inference use a fixed internal seed (`12345`).
