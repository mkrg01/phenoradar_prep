# Updating older configurations and results

[Documentation](index.md) · [Current configuration](configuration.md)

Use this guide only when resuming an analysis created with an older workflow.
Compare your overrides with [config/config.yaml](../config/config.yaml). The
workflow rejects the obsolete settings below rather than silently ignoring them.

## Configuration changes

Rename `analysis` to `run_name` (default `run001`). For default container
execution, change `container_image: null` to `auto` or a matching image path/URI.
For native execution, pass `--software-deployment-method conda`.

Generated paths and tool commands are now fixed. Remove these old overrides:

| Section | Keys to remove |
| --- | --- |
| Top level | `tools`, `paths` |
| `taxonomy` | `database` |
| `odb` | `version`, `reference_dir` |
| `kegg` | `command`, `reference_dir` |
| `alignment` | `command`, `famsa_command` |
| `phylogeny` | `famsa_command`, `trimal_command`, `veryfasttree_command`, `astral_command` |
| `phylogeny.dating` | `command`, `threads` |
| `phylogeny.dating.timetree` | `python`, `request_delay_seconds`, `cache_dir` |

Other changes:

- Replace `odb.mem_mb` with `odb.mem_gb`, dividing by 1000. Cache storage is now
  `.cache/`; `PHENORADAR_CACHE_DIR` is unused.
- Marker selection ranks overall coverage and caps loci. Remove `phylogeny` keys
  `clade_rank`, `clade_min_species`, `min_clade_occupancy`, `min_occupancy`,
  `max_gap_fraction`, `min_sites`, `min_informative_sites`, and `min_markers_per_species`.
  Use `max_markers`, `trimal_mode`, and per-sequence `min_protein_length`.
- Replace `phylogeny.dating.treepl` with `lsd2`; rate partitions are unsupported.
  Remove `phylogeny.dating.timetree.min_clade_taxa`.
- Taxonomy review uses MonoPhy. Remove `taxonomy_audit.min_reference_species`
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
./run_pipeline.sh --configfile config/mydata.yaml \
  --cores 1 --resources mem_gb=4 -- phenoradar_metadata
```

This uses completed selected metadata and the trait source. See
[PhenoRadar inputs](phenoradar_inputs.md) for the collection requirements.

## Relocated results

Update absolute CDS/abundance paths in sample manifests after moving inputs.
Refresh PhenoRadar links and regenerate filtered exports from relocated sources.
Verify the new paths, preserve historical provenance, and record the relocation
separately. See [outputs](outputs.md) for the current layout.
