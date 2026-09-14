# Trait contrast pairs

[Documentation](index.md) · [BUSCO phylogeny](phylogeny.md)

Pairs use nwkit's homogeneous-clade grouping and contrastive-clade selection.
Traits come from `inputs.species_trait`, selected by `contrast.trait` (default `C4`).
Missing traits remain unknown; multiple expression runs do not add species weight.

| Target | Tree | Species membership |
| --- | --- | --- |
| `phylogeny_contrast_pairs` | Full/phenotyped BUSCO tree | Assigned directly on the molecular tree |
| `contrast_pairs` | Newly inferred representative tree | Inherited through NCBI and molecular group maps |

## Pairs from full or phenotyped trees

Choose `phylogeny.species_sets: [all]`, `[phenotyped]`, or both:

```bash
./run_pipeline.sh --configfile config/mydata.yaml \
  --cores 32 --resources mem_gb=128 -- phylogeny_contrast_pairs
```

Without exclusions, missing or outdated trees can schedule inference. Pair-only
jobs need one core and 8 GB; confirm with `--dry-run` before reducing the budget.
`contrast.enabled` does not control this target.

Inputs are the tree/QC, sample manifest, selected BUSCO metadata, and traits.
Missing-trait tips are pruned while preserving root direction and path lengths;
`observed_tree.nwk` records the subtree. Unknown-trait species remain in output
metadata with empty pair IDs. Zero/one observed state gives zero pairs; more than
two states is unsupported.

Results are in `phylogeny/all/contrast/` or `phylogeny/phenotyped/contrast/` under
`results/<run_name>/`.

## Representative analysis

```bash
./run_pipeline.sh --configfile config/mydata.yaml \
  --cores 32 --resources mem_gb=128 -- contrast_pairs
```

This requires the ordinary [phylogeny inputs](phylogeny.md#inputs), exactly two
observed states, and at least four representatives. Set `contrast.enabled: true`
to include it in `all`; outputs go to `phylogeny/representatives/`.

1. Skim homogeneous clades on the local NCBI guide, choosing representatives by
   complete BUSCO fraction; `phylogeny.seed` controls randomized ties.
2. Resolve the outgroup within that subset and infer its BUSCO tree independently.
3. Skim the molecular tree and pair minimal mixed clades with two opposite-state
   representatives. Multiway cases remain unresolved.
4. Map original species through both grouping stages to final representatives/pairs.

A manual outgroup must be selected among the representatives. Dating is not part
of this target. To inspect selection alone, target
`results/<run_name>/phylogeny/representatives/selection/selection.json`.

## After species exclusion

With nonempty `exclude_species`, `phylogeny_contrast_pairs` uses the
[filtered export](species_filter.md) and requires completed inputs for the
requested sets. It never schedules inference in this mode.

Filtering recomputes pairs for completed full and phenotyped trees, regardless
of `phylogeny.species_sets`, using their original root direction. Results go to
`filtered/phylogeny/<set>/contrast/`. Representatives are not filtered.

## Outputs

Each branch's `contrast/` directory contains:

| Output | Contents |
| --- | --- |
| `contrast_pairs.tsv` | Pair states, representatives, and species counts |
| `species_metadata.tsv` | Species traits, groups, representatives, and nullable pair IDs |
| `summary_tree.nwk`, `.all.tsv`, `.sampled.tsv` | Molecular skim and membership |
| `contrastive.nwk`, `.all.tsv`, `.sampled.tsv` | Raw contrastive selection, including unresolved candidates |
| `summary_tree.pdf`, `summary_tree.svg` | Trait-colored tree with group/pair IDs |
| `summary.json` | Counts, settings, and provenance |

Representative selection records are in `selection/`. Zero-pair results are valid
and have header-only pair tables; empty trees are recorded when no tips qualify.
To restore a deleted file inside a filtered bundle, run `filter_species` with
`--forcerun filter_species`.

## Interpretation

Pairs describe sampled trait contrasts, not independent evolutionary origins.
NCBI-inherited membership has not been molecularly tested for every species.
Assignments are not filtered by branch support. Pair IDs are local to each result
and can change with species membership or exclusions.
