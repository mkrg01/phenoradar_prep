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
| `summary_tree.pdf`, `.svg` | Trait-colored tree with group/pair IDs |
| `summary.json` | Counts, settings, and provenance |

Pairs describe sampled trait contrasts, not independent evolutionary origins.
NCBI group membership is not molecularly tested for every species, and selection
uses no branch-support filter. `n_species_*` counts biological species;
`n_samples_*` counts sample tips.

Pair IDs are local to each result and may change with sample selection. Normal
analysis exclusions apply before inference. [Post hoc filtering](species_filter.md)
can recompute pairs for completed all/phenotyped trees, but cannot export
representative results. Zero-pair results are valid.
