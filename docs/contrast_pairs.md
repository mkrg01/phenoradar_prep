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
  trees: [all]
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
| `expanded_pair_clades.tsv` | Anchoring pairs, expanded-clade membership counts, and both trait-state counts |
| `summary_tree.nwk`, `.all.tsv`, `.sampled.tsv` | Grouped molecular tree and membership |
| `contrastive.nwk`, `.all.tsv`, `.sampled.tsv` | Raw selection, including unresolved candidates |
| `summary_tree.pdf`, `.svg` | Publication view: scientific names, group species counts, and red pair clades |
| `summary_tree_id.pdf`, `.svg` | Same tree with run, group, and pair IDs in tip labels |
| `expanded_pair_tree.pdf`, `.svg` | Same tree with expanded-clade shading and blue expanded branches; original pair branches remain red |
| `expanded_pair_tree_id.pdf`, `.svg` | Expanded view with run, group, original pair, and expanded-clade IDs |
| `summary.json` | Counts, settings, and provenance |

For PhenoRadar input, use `metadata/species_metadata.tsv` at the analysis root.
The `all` target joins these pair assignments to the cohort's traits and taxonomy.
When several tree branches are selected, each has a separate pair-ID column in
that table. See [input metadata](phenoradar_inputs.md#input-files).

Figures show the trait states (`0` and `1`) and biological species counts per group.
Use the ID views to trace representatives, groups, and pairs back to sample
records. Red branches mark original pairs; the expanded views also mark the
larger clades.

Pair IDs describe sampled contrasts. Inferring the number of evolutionary
origins requires a separate analysis. Pair selection uses topology and trait
states without a branch-support cutoff. On the representative tree, members
inherit their representative's placement from the NCBI grouping. Counts
distinguish biological species (`n_species_*`) from sample tips (`n_samples_*`).

[Filtering a completed tree](species_filter.md) can recompute pairs for
all/phenotyped results. Pair IDs are specific to each tree and cohort.

## Expanded pair clades

`expanded_pair_clade_id` assigns the largest rooted clade containing one selected
contrast pair and no other selected pair. Start at each pair's common ancestor
on the molecular `summary_tree.nwk`, then ascend until the parent contains
another pair. Unpaired tips remain available during expansion; do not use the
pair-only `contrastive.nwk` for this step.

Expanded clades are disjoint and use the anchoring pair IDs. An added sample
can have an expanded-clade ID and a blank original pair ID. Missing-trait
samples and branches outside the expanded clades have blank IDs. With a single
pair, its expanded clade can include the entire observed tree.

For leave-one-group-out validation, review the number of clades and their
sample counts and label balance. An expanded clade can be much larger than
its anchor pair. See [PhenoRadar validation groups](phenoradar_inputs.md#validation-groups).

### Apply to saved results

Existing completed analyses are immutable. Export a derivative from saved pairs
without inference or pair reselection:

```bash
python workflow/scripts/expanded_pair_clades.py --analysis results/leaf/downstream/carnivory
```

Run with the workflow test environment or the pinned `timetree` environment.
Every configured contrast branch must already have its summary tree, metadata,
pair table, and summary. Outputs go under `<analysis>/expanded_pair_clades/`:

- `metadata/species_metadata.tsv`: the full cohort, taxonomy, original pairs,
  and expanded-clade columns; its JSON records the joins.
- `phylogeny/<tree>/contrast/`: the unchanged saved tree and pair table,
  augmented branch metadata, expanded-clade counts, and both expanded figures.
- `provenance.json`: source/output checksums, code version, and branch counts.

Expression inputs remain at their original analysis paths. Publication is
atomic and existing export directories are refused; use `--outdir NEW_DIRECTORY`
for another export. Original files and completion records are unchanged.
