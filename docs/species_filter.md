# Manual species exclusion

[Documentation](index.md) · [PhenoRadar inputs](phenoradar_inputs.md)

`filter_species` exports completed results to `results/<run_name>/filtered/`,
removing all runs and gene copies of named species while preserving the original
analysis. Choose exact IDs from `metadata/samples.tsv` in a top-level YAML list:

```yaml
exclude_species:
  - Lespedeza_davurica
  - Cleistogenes_squarrosa
```

These are format examples. The default is `[]`; taxonomy check flags never set
exclusions automatically. Unknown/duplicate IDs and removal of all species fail.

## Execute after the source analysis

Save exclusions in the dataset config or a local override:

```bash
./run_pipeline.sh --configfile config/mydata.yaml config/exclusions.local.yaml \
  --cores 1 --resources mem_gb=8 -- filter_species
```

The target needs the original sample manifest and completed outputs, not raw
sequences or reference preparation. It exports complete branches and reports
absent/incomplete branches in `manifest.json`; it never starts upstream analyses.
Inconsistent identities/checksums or missing files claimed by completion records fail.

Changing exclusions or source outputs refreshes the export. Every export starts
from the original analysis, so removing an exclusion restores that species.
For archived results outside the workflow layout:

```bash
python workflow/scripts/filter_species.py \
  --source /path/to/completed/analysis \
  --exclude-species '["Lespedeza_davurica"]' \
  --traits input/species_trait.tsv --outdir /path/to/new/filtered
```

`--outdir` defaults to `<source>/filtered`. The script requires the workflow's
`timetree` environment; it makes no TimeTree query.

## Exported dataset

Only available, completed sections are exported:

| Output under `filtered/` | Contents |
| --- | --- |
| `metadata/`, `excluded_samples.tsv` | Retained metadata and excluded run/species identities |
| `proteins/` | Symlinks to retained species' proteins |
| `orthogroups/mapping/` | Retained gene ownership and all gene/OG assignments |
| `orthogroups/expression/`, `kegg/` | Expression, membership, support, and QC with excluded runs/genes removed |
| `orthogroups/alignments/` | Retained gene rows; empty OGs omitted and listed in `filter_qc.json` |
| `phylogeny/all/` | Pruned species/gene trees, available alignments, and dated tree |
| `phylogeny/all/contrast/`, `phylogeny/phenotyped/contrast/` | Recomputed pairs when tree/QC, manifest, BUSCO scores, and traits are ready |
| `manifest.json` | Available/skipped sections, exclusions, identities, and checksums |

OG alignments require a completed inventory and valid [gene IDs](alignments.md#outputs-and-phenoradar).
The representative analysis, phenotyped inference files, and taxonomy check reports
stay in the source analysis. Completed phenotyped trees can still supply new pairs.
Missing pair inputs are recorded, without triggering inference.

Pair assignment uses `contrast.trait` and top-level `seed` (standalone options
`--contrast-trait` and `--seed`). Exclusions can change clades and pair IDs;
zero/one-state subsets produce zero pairs. See [contrast pairs](contrast_pairs.md).

## Numerical and phylogenetic meaning

Expression values and feature axes are retained without reaggregation or
normalization. KO zeros remain distinct from unavailable values. Alignments keep
all columns, including newly all-gap columns; marker selection and alignment
are not repeated.

Pruned trees retain path lengths but are not new inference or dating estimates.
Internal supports are omitted because they have not been recalculated. If the
outgroup is removed, `pruning.json` marks the root for review; no new root is
chosen. Original node-age/calibration tables remain with the original analysis.

## Storage and verification

Tables and databases take new disk space. Protein symlinks and metadata paths
still refer to original files; keep those sources available. Use the
[collector](phenoradar_inputs.md#species-exclusions) to pass the subset to PhenoRadar.
