# Trait contrast pairs

[Documentation](index.md) · [BUSCO phylogeny](phylogeny.md)

Select opposite-trait representatives with nwkit's homogeneous-clade grouping
and contrastive-clade selection. Traits come from the column named by `trait`
in `inputs.species_trait`.

## Configuration and execution

Before preparing the analysis:

```yaml
trait: carnivory
phylogeny:
  trees: [representatives]
  contrast_pairs:
    enabled: true
```

```bash
./run_analysis.sh submit --analysis results/leaf/downstream/carnivory --target contrast_pairs
```

The target includes missing tree inference and requires the enabled flag.
Pairs are also included in `all`.

| Tree | Selection before pairing |
| --- | --- |
| `all` | Infer all selected samples, then prune missing-trait tips |
| `phenotyped` | Remove missing-trait samples before inference |
| `representatives` | Remove missing-trait samples and select NCBI representatives before inference |

All/phenotyped trees with zero or one observed state produce zero pairs.
More than two states is unsupported. Missing-trait samples keep empty pair IDs.

## Representative selection

Representative inference requires exactly two observed states and at least four
representatives. It also works with pair selection disabled.

1. Group homogeneous trait clades on the NCBI guide and choose representatives
   by BUSCO completeness, breaking ties with seed `12345`.
2. Resolve the outgroup and infer a BUSCO tree for those representatives.
3. Group the molecular tree and pair minimal mixed clades with opposite-state
   representatives; unresolved multiway cases remain unpaired.
4. Map original samples through the groups to their representatives and pairs.

Both trait states are reduced. A manual outgroup must belong to the selected
representatives. Review `phylogeny/representatives/selection/`.

## Outputs and interpretation

Under `results/<build>/downstream/<analysis>/phylogeny/<tree>/contrast/`:

| Output | Contents |
| --- | --- |
| `contrast_pairs.tsv` | Pair states, representatives, and counts |
| `species_metadata.tsv` | Traits, groups, representatives, and nullable pair IDs |
| `summary_tree.nwk`, `.all.tsv`, `.sampled.tsv` | Grouped molecular tree and membership |
| `contrastive.nwk`, `.all.tsv`, `.sampled.tsv` | Raw selection, including unresolved candidates |
| `summary_tree.pdf`, `.svg` | Publication view: scientific names, group species counts, and red pair clades |
| `summary_tree_id.pdf`, `.svg` | Same tree with run, group, and pair IDs in tip labels |
| `summary.json` | Counts, settings, and provenance |

For PhenoRadar input, use `metadata/species_metadata.tsv` at the analysis root.
The `all` target joins these pair assignments to the cohort's traits and taxonomy.
When several tree branches are selected, each has a separate pair-ID column in
that table. See [input metadata](phenoradar_inputs.md#input-files).

The legend shows the configured trait name and source state labels. With
[explicit text encoding](inputs.md#traits), the label mapped to `0` is gray
and the label mapped to `1` is green; numeric inputs display `0`/`1`, including
for `C4`. `summary.json` records `trait_labels` and `state_order`, preserving the
labels and colors in post hoc subsets. Older summaries without labels use their
state values. The plot reserves a separate legend area and expands its width
for long names.

Tip labels show the representative's scientific name in italics, followed by
upright `(n = ...)`, where `n` counts biological species in its group. Red branches
mark the two-tip clade of each selected pair; unpaired branches remain black.
The ID view adds run IDs from `metadata/samples.tsv`, plus group and pair IDs
(`-` means unpaired). It retains sample IDs when explicit run information is
unavailable. Only scientific names are italicized. Both views use the same
topology, branch lengths, trait colors, and pair assignments; drawing requires
no inference or pair reselection. Both views remain in the contrast output directory.

Pairs describe sampled trait contrasts, not independent evolutionary origins.
NCBI group membership is not molecularly tested for every species, and selection
uses no branch-support filter. `n_species_*` counts biological species;
`n_samples_*` counts sample tips.

Pair IDs are local to each result and may change with sample selection. Normal
analysis exclusions apply before inference. [Post hoc filtering](species_filter.md)
can recompute pairs for completed all/phenotyped trees, but cannot export
representative results. Zero-pair results are valid.
