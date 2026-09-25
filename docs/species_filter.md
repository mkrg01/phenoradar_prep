# Manual species exclusion

[Documentation](index.md) · [PhenoRadar inputs](phenoradar_inputs.md)

Normally, set `exclude_species` in `config/analysis.yaml` **before preparing a new
analysis**. Those species are removed before any computation, including tree
inference. Use exact species IDs from the build; unknown/duplicate IDs or an empty
selection are errors. Taxonomy flags never exclude species automatically.

```yaml
exclude_species:
  - Lespedeza_davurica
  - Cleistogenes_squarrosa
```

To exclude an RNA-seq run from future builds as well, use the
[accession list](datasets.md#manually-excluding-unusable-accessions).

## Execute after the source analysis

Alternatively, `filter_species` exports a subset of completed results without
repeating analysis. Save the YAML above as `config/export_exclusions.yaml`, using
IDs from the source `metadata/samples.tsv`, and invoke the low-level launcher:

```bash
./run_pipeline.sh --cores 1 --resources mem_gb=8 \
  --configfile analyses/analysis001/pipeline.yaml config/export_exclusions.yaml -- filter_species
./run_pipeline.sh --cores 1 --resources mem_gb=4 \
  --configfile analyses/analysis001/pipeline.yaml config/export_exclusions.yaml -- phenoradar_inputs
```

Use the same override for both commands; do not edit the frozen `pipeline.yaml`.
Filtering writes `results/<analysis>/filtered/`, preserves the original analysis,
and never starts producer jobs. It validates completed branches and reports
missing/incomplete ones in `manifest.json`. Rerunning refreshes from the originals;
removing an exclusion restores that species. To recreate deleted export files,
add `--forcerun filter_species` before `--`.

## Exported dataset

| Output under `filtered/` | Contents |
| --- | --- |
| `metadata/`, `excluded_samples.tsv` | Retained metadata and excluded identities |
| `proteins/` | Links to retained proteins |
| `orthogroups/`, `kegg/` | Filtered mappings, expression, alignments, support, and QC |
| `phylogeny/all/` | Pruned species/gene trees, alignments, and dated tree |
| `phylogeny/{all,phenotyped}/contrast/` | Recomputed pairs when required source records are complete |
| `manifest.json` | Exported/skipped sections, exclusions, and checksums |

Representative results, phenotyped inference files, and taxonomy reports remain
in the source analysis. Completed phenotyped trees can still supply new pairs.
Pairs use `trait` and `seed`; pair IDs can change after exclusion.

## Numerical and phylogenetic meaning

Expression values and feature axes are preserved without reaggregation or
normalization, including KO zeros versus unavailable values. Alignments retain
all columns; empty OGs are omitted and listed in `filter_qc.json`.

Pruned trees retain path lengths but are not new inference/dating estimates.
Internal supports are removed. Removing the outgroup marks the root for review
in `pruning.json`; no new root is chosen. Original age/calibration tables remain
with the source. Keep source files available: protein links and metadata paths
depend on them.
