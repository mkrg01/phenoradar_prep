# Updating older configurations and results

[Documentation](index.md) · [Current configuration](configuration.md)

Use this guide only when resuming an analysis created with an older workflow.
Compare your overrides with [config/config.yaml](../config/config.yaml). The
workflow rejects the obsolete settings below rather than silently ignoring them.

## Configuration changes

Rename the top-level `analysis` key to **`run_name`**. The default is now
`run_name: run001`, which selects `results/run001/`, `work/run001/`, and
`logs/run001/`. Command-line overrides use `--config run_name=YOUR_RUN_NAME`.

`run_pipeline.sh` now enables Singularity in both direct and Slurm execution.
The default `container_image: auto` selects the image matching `VERSION`.
Existing dataset configurations with `container_image: null` must change it to
`auto` to use this default, or specify a matching image URI/absolute SIF path;
see [container setup](containers.md). For native Conda execution, pass
`--software-deployment-method conda`; either `auto` or `null` works.

Generated paths, tool commands, OrthoDB v12, serial LSD2 execution, and the
one-second delay between uncached TimeTree requests are fixed by the workflow.
Remove these old override keys:

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

Replace `odb.mem_mb` with `odb.mem_gb`, dividing by 1000: `128000` MB becomes
`128` GB. The launcher now keeps its cache in `.cache/`;
`PHENORADAR_CACHE_DIR` is no longer used. `SNAKEMAKE_BIN` remains available for
selecting the launcher executable.

Phylogeny marker selection now ranks by overall coverage and caps the number
of loci. Remove `phylogeny.clade_rank`, `clade_min_species`,
`min_clade_occupancy`, and `min_occupancy`. Remove `max_gap_fraction`; column
selection uses `trimal_mode`. The former `min_sites`, `min_informative_sites`,
and `min_markers_per_species` settings are also removed. Counts remain in QC,
while the per-sequence `min_protein_length` filter still applies.

Remove `phylogeny.dating.timetree.min_clade_taxa`: small clades are eligible,
with workload bounded by representative and query caps. Replace
`phylogeny.dating.treepl` with the current `lsd2` settings. Rate partitions are
unsupported. Old treePL CV/restart files are historical outputs; current LSD2
provenance identifies the completed dating result.

Taxonomy review uses MonoPhy. Remove `taxonomy_audit.min_reference_species`
and `max_plot_species`; the old custom detector's envelopes, gene-tree scores,
and context plots are no longer generated. Regenerating an owned audit directory
replaces the report, including obsolete files.

Use a top-level `exclude_species` list rather than a `species_filter` wrapper.
Exclusions apply to every run of each named species.

The KEGG ambiguity default is `duplicate`. An explicit older `drop` override
still drops multi-KO genes; change it only if full TPM contributions to every
accepted KO are intended. See [KO quantification](kegg.md#assignment-and-quantification).

## Importing ODB results

`odb.existing_results` imports a verified snapshot of earlier annotations:

```yaml
odb:
  existing_results: /path/to/odb_snapshot
```

The directory must contain `annotations.tsv` and `snapshot.json` with schema
version 1, OrthoDB version/node, original protein checksums keyed by species,
and the annotation checksum. Snapshot preparation is a local migration step;
there is no workflow target that reconstructs evidence for an old run.
The validation contract is in [merge_odb.py](../workflow/scripts/merge_odb.py).

With this setting, `mapping` and `all` skip ODB reference preparation and mapping.
Translation still verifies exact agreement with the original protein inputs.
Missing species, changed proteins, mismatched version/node, or modified
annotations stop the import. Species subsets are supported; rows outside the
selected FASTA IDs are removed and counted in `merge_qc.json`. The normal
merged mapping and TPM outputs are produced. The explicit `references` target
still prepares OrthoDB when requested.

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

After moving inputs, update absolute CDS and abundance paths in selected-sample
manifests. Refresh PhenoRadar links after moving producer outputs, and regenerate
filtered exports from the relocated sources before reuse. Verify files and links
at their new locations.

Preserve historical `run.json` and provenance records. Archive the original
operational files and record path changes separately; the earlier local migration
used `logs/layout_migration/<timestamp>/`. See [outputs](outputs.md) for the current layout.
