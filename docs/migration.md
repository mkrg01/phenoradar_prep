# Updating older configurations and results

[Documentation](index.md) · [Current configuration](configuration.md)

Use this guide only when resuming an analysis created with an older workflow.
Compare your overrides with [config/config.yaml](../config/config.yaml). The
workflow rejects the obsolete settings below rather than silently ignoring them.

## Configuration changes

Rename `analysis` to `run_name` (default `run001`). Change `container_image` values
of `null`, image URIs, or SIF paths to `auto` or a GHCR release version such as
`"0.1.0"` (without `v`). Only GHCR release versions are supported.
For native execution, pass `--software-deployment-method conda`.

Generated paths and tool commands are now fixed. Remove these old overrides:

| Section | Keys to remove |
| --- | --- |
| Top level | `tools`, `paths`, `taxonomy` (including `source` and `database`) |
| `odb` | `version`, `reference_dir`, `min_free_gb`, `reference_min_free_gb`, `allow_nonlocal`, `keep_work` |
| `kegg` | `command`, `reference_dir` |
| `alignment` | `command`, `famsa_command` |
| `phylogeny` | `famsa_command`, `trimal_command`, `veryfasttree_command`, `astral_command`, `sequence_mode`, `sequence_dir`, `sequence_suffix`, `busco_full_suffix` |
| `phylogeny.dating` | `command`, `threads` |

Other changes:

- Move `phylogeny.seed` to top-level `seed`, keeping its value. The default is
  still `12345`; the seed is shared by tree inference and representative/pair
  selection. The old nested key is rejected. See [reproducibility](configuration.md#reproducibility).
- Rename `taxonomy_audit` to `taxonomy_check` in configuration and targets.
  The report rule is now `check_taxonomy`; reports and logs use `taxonomy_check`.
  Run the new target to create reports at the new paths, reusing completed
  species trees. Earlier report directories are left intact.
- Remove the entire `phenoradar` section. `phenoradar_inputs` now automatically
  collects available completed outputs, including both KEGG grouping maps and
  all tree/contrast branches. Trees and contrast metadata keep their original
  branch paths; no single tree or pair assignment is selected. Collection also
  accepts multiple runs per species and datasets without expression results.
  See [collected inputs](phenoradar_inputs.md) for the layout.
- Remove all CPU and memory keys from workflow config, including `threads`,
  `mem_gb`, and phylogeny's `*_threads` / `*_mem_gb`. Defaults now live in each
  rule. Use `--set-threads` and `--set-resources` for
  [per-rule overrides](running.md#resource-budgets).
- Remove `odb.chunk_size` and `odb.batch_size`, including from pilot overrides.
  The rules now use 100 species per chunk and an internal batch size of four
  times the allocated ODB threads. Existing smaller chunks will be regrouped;
  mapping work is reused only when chunk contents and execution settings match.
- Phylogeny always uses the selected samples' CDS from `inputs.cds_dir`. Move any
  previous CDS `sequence_dir` override there. BUSCO full-table filenames and gzip
  versions are detected automatically; use the [supported layouts](phylogeny.md#inputs).
- Taxonomy downloads automatically when `resources/taxonomy/taxa.sqlite` is
  missing. Existing snapshots are reused across runs; see [references](references.md#taxonomy-reference).
- Storage capacity and filesystem choice are left to the user. ODB and KofamScan
  work files are retained after success; see [recovery](running.md#re-running-and-recovery).
- `phylogeny.dating.calibration_source` defaults to `timetree`. Set it to `file`
  to keep using your own calibration TSV.
- Remove the entire `phylogeny.dating.lsd2` section. Dating now always uses
  input-length weights, an automatic variance offset, and the sum of retained
  trimmed alignment lengths as the site count. The standalone helper's
  `--settings` option is also removed; fixed choices remain recorded in dating
  provenance. See [LSD2 fitting](dating.md#lsd2-fitting).
- Remove the entire `phylogeny.dating.timetree` section, including
  `max_representatives`, `max_queries`, and `min_studies`. Every internal node
  is now considered, with all resolved descendant species and a fixed minimum
  of five distinct bibliographic records. The helper no longer accepts
  `--settings`, `--coverage`, or `--representatives`. Cached responses remain
  compatible and are reused; only missing responses are fetched. Regenerating
  calibrations replaces `representatives.txt` / `representatives.nwk` with
  `nodes.nwk` and adds `studies.tsv`. Completed species-tree inference is reused.
  See [TimeTree calibrations](dating.md#timetree-calibrations).
- `selection.species_list` now specifies candidates before BUSCO filtering;
  listed species below the threshold or without BUSCO data are excluded.
- Remove the older `odb.mem_mb` override too. Cache storage is now `.cache/`;
  `PHENORADAR_CACHE_DIR` is unused.
- Marker selection ranks overall coverage and caps loci. Remove `phylogeny` keys
  `clade_rank`, `clade_min_species`, `min_clade_occupancy`, `min_occupancy`,
  `max_gap_fraction`, `min_sites`, `min_informative_sites`, and `min_markers_per_species`.
  Use `max_markers`, `trimal_mode`, and per-sequence `min_protein_length`.
- Remove `phylogeny.dating.treepl`; dating uses LSD2 with no rate partitions.
- Taxonomy review uses MonoPhy. Remove `taxonomy_check.min_reference_species`
  and `max_plot_species`. Regeneration replaces the old report.
- Replace the `species_filter` wrapper with top-level `exclude_species`.
- KEGG ambiguity now defaults to `duplicate`. Explicit `drop` remains valid;
  see [KO quantification](kegg.md#assignment-and-quantification) before changing it.

## Importing ODB results

`odb.existing_results` imports a verified snapshot of earlier annotations:

```yaml
odb:
  existing_results: /path/to/odb_snapshot
```

The snapshot requires `annotations.tsv` and schema-v1 `snapshot.json` recording
OrthoDB version/node, per-species original protein checksums, and the annotation
checksum. See [merge_odb.py](../workflow/scripts/merge_odb.py) for the format.
There is no target to reconstruct missing provenance for an old run.

`mapping`/`all` reuse annotations and skip ODB preparation/mapping. Translation
still verifies proteins; missing species, changed proteins, and mismatched
reference/annotation checksums fail. Subsets are supported. The explicit
`references` target still prepares OrthoDB.

## Backfilling PhenoRadar metadata

If completed results lack `metadata/species_metadata.tsv`, prepare it before
collecting downstream inputs:

```bash
./run_pipeline.sh --cores 1 --resources mem_gb=4 -- phenoradar_metadata
```

This uses completed selected metadata and the trait source. See
[PhenoRadar inputs](phenoradar_inputs.md) for the collection requirements.

## Relocated results

Update absolute CDS/abundance paths in sample manifests after moving inputs.
Refresh PhenoRadar links and regenerate filtered exports from relocated sources.
Verify the new paths, preserve historical provenance, and record the relocation
separately. See [outputs](outputs.md) for the current layout.
