# Exclude species or export a subset

[Documentation](index.md) · [PhenoRadar inputs](phenoradar_inputs.md)

For a new analysis, set `exclude_species` in `config/analysis.yaml`. Biological
species IDs exclude all their samples; exact sample IDs exclude only that sample.
Exclusions apply before computation, including tree inference.

```yaml
exclude_species:
  - Lespedeza_davurica
  - Cleistogenes_squarrosa
```

To exclude runs from future builds, use the
[accession list](datasets.md#manually-excluding-unusable-accessions).

## Export from completed results

`filter_species` creates a subset without repeating the analysis. Save an
`exclude_species` list in `config/export_exclusions.yaml`, using exact sample IDs
from the source `metadata/samples.tsv`. This low-level command does not expand
biological species IDs; list every sample to exclude.

```bash
./run_pipeline.sh --cores 1 --resources mem_gb=8 \
  --configfile results/leaf/downstream/carnivory/pipeline.yaml config/export_exclusions.yaml -- filter_species
./run_pipeline.sh --cores 1 --resources mem_gb=4 \
  --configfile results/leaf/downstream/carnivory/pipeline.yaml config/export_exclusions.yaml -- phenoradar_inputs
```

Use the same override for both commands. The export goes to `filtered/` under
the source analysis, which stays unchanged. No producer jobs are started;
`manifest.json` lists exported and incomplete sections. Rerunning refreshes from
the original results. To recreate deleted exports, add `--forcerun filter_species`
before `--`.

## What the export contains

| Output under `filtered/` | Contents |
| --- | --- |
| `metadata/`, `excluded_samples.tsv` | Retained metadata and excluded samples |
| `proteins/` | Links to retained proteins |
| `orthogroups/`, `kegg/` | Subset mappings, expression, alignments, and QC |
| `phylogeny/all/` | Pruned species/gene trees, alignments, and dated tree |
| `phylogeny/{all,phenotyped}/contrast/` | Recomputed pairs when source records are complete |
| `manifest.json` | Exported/skipped sections and checksums |

Representative results, phenotyped inference files, and taxonomy reports stay
in the source. Completed phenotyped trees can still supply new pairs.

Expression values and feature axes are preserved, including KO zeros versus
unavailable values. KO labels retain their original representative metadata,
even if a donor sample is excluded. Alignments keep all columns; empty OGs are
omitted and listed in `filter_qc.json`.

Pruned trees preserve path lengths but are not new inference or dating estimates.
Internal supports are removed. Removing the outgroup marks the root for review
in `pruning.json`; no new root is chosen. Pair IDs may change. Keep source files
available for linked proteins, metadata, and original calibration records.
