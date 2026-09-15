# Manual species exclusion

[Documentation](index.md) · [PhenoRadar inputs](phenoradar_inputs.md)

`filter_species` exports completed results to `results/<run_name>/filtered/`,
removing all runs and gene copies of named species while preserving the original
analysis. Set exact IDs from `metadata/samples.tsv` in `config/config.yaml`, e.g.:

```yaml
exclude_species:
  - Lespedeza_davurica
  - Cleistogenes_squarrosa
```

The default is `[]`. Unknown/duplicate IDs or removing all species are errors;
taxonomy flags never set exclusions automatically.

## Execute after the source analysis

```bash
./run_pipeline.sh --cores 1 --resources mem_gb=8 -- filter_species
```

This requires the original sample manifest and completed outputs. It exports
complete branches and reports absent/incomplete ones in `manifest.json`, without
starting analyses. Invalid identities, checksums, or missing recorded files fail.

Rerunning refreshes changed exclusions/results from the original analysis;
removing an exclusion restores that species. To restore deleted export files,
use `--forcerun filter_species`. For archived analyses, see
[filter_species.py](../workflow/scripts/filter_species.py) `--help` in the
workflow's `timetree` environment.

## Exported dataset

| Output under `filtered/` | Contents |
| --- | --- |
| `metadata/`, `excluded_samples.tsv` | Retained metadata and excluded identities |
| `proteins/` | Symlinks to retained species' proteins |
| `orthogroups/`, `kegg/` | Filtered mappings, expression, alignments, support, and QC |
| `phylogeny/all/` | Pruned species/gene trees, available alignments, and dated tree |
| `phylogeny/{all,phenotyped}/contrast/` | Recomputed pairs when tree/QC, manifest, BUSCO scores, and traits are complete |
| `manifest.json` | Exported/skipped sections, exclusions, and checksums |

Representative analysis, phenotyped inference files, and taxonomy reports remain
in the source analysis. Completed phenotyped trees can still supply new pairs.
Pairs use `contrast.trait` and top-level `seed`; IDs can change after exclusion,
and zero/one-state subsets yield zero pairs.

## Numerical and phylogenetic meaning

Expression values and feature axes are preserved without reaggregation or
normalization, including KO zeros versus unavailable values. Alignments retain
all columns; empty OGs are omitted and listed in `filter_qc.json`.

Pruned trees retain path lengths but are not new inference/dating estimates.
Internal supports are removed. Removing the outgroup marks the root for review
in `pruning.json`; no new root is chosen. Original age/calibration tables stay
with the source analysis.

Tables/databases use new disk space; protein symlinks and metadata paths depend
on original files. Keep sources available and use the
[collector](phenoradar_inputs.md#species-exclusions) to pass the subset downstream.
