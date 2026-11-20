# Contrasting sister-clade pairs and phylogenetic blocks

[Documentation](index.md) · [BUSCO phylogeny](phylogeny.md)

Select pairs of sister clades with different states of the focal trait using
nwkit. Traits come from the column named by `trait` in `inputs.species_trait`.

- A **trait clade** is a maximal clade with one observed trait state among its
  sampled members; it can contain one or several species.
- A **contrast pair** consists of two sister trait clades with opposite states.
  All members of both clades retain pair membership.
- A **phylogenetic block** is the largest rooted clade containing one selected
  contrast pair and no other selected pair. It may include unpaired samples.

Public metadata uses `trait_clade_id`, `contrast_pair_id`, and
`phylogenetic_block_id`. The raw nwkit skim tables retain nwkit's `group` column.
Readers accept the former `group` and `expanded_pair_clade_id` metadata columns
and the former `expanded_pair_clades.tsv` filename. Conflicting old and new
values are errors. New outputs use the canonical names.

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

For `all` and `phenotyped`, pairs are identified directly on the molecular tree
after removing missing-trait tips. Both sister clades must be uniform for the
trait and have opposite states. Blocks are also assigned on that tree, before
selecting any display representatives. BUSCO completeness and the random seed
affect display representatives, not pair IDs or pair/block membership.

## Representative selection

Representative inference requires exactly two observed states and at least four
representatives. It also works with pair selection disabled.

1. Group homogeneous trait clades on the NCBI guide and choose representatives
   by BUSCO completeness, breaking ties with seed `12345`.
2. Resolve the outgroup and infer a BUSCO tree for those representatives.
3. Identify opposite-state sister trait clades and their blocks directly on the
   inferred molecular tree; unresolved multiway cases remain unpaired.
4. Map original samples through the groups to their representatives and pairs.

Both trait states are reduced. A manual outgroup must belong to the selected
representatives. Review `phylogeny/representatives/selection/`.

After assigning pairs and blocks, summary figures retain one sample per trait
clade, selected by BUSCO completeness with seeded tie breaking. These display
samples do not change the assignments. The `representatives` inference route
still requires its separate NCBI-based selection before tree estimation.

## Outputs and interpretation

Under `results/<build>/downstream/<analysis>/phylogeny/<tree>/contrast/`:

| Output | Contents |
| --- | --- |
| `contrast_pairs.tsv` | Pair states, representatives, and counts |
| `species_metadata.tsv` | Traits, groups, representatives, and nullable pair IDs |
| `phylogenetic_blocks.tsv` | Anchoring pairs, expanded-clade membership counts, and both trait-state counts |
| `summary_tree.nwk`, `.all.tsv`, `.sampled.tsv` | Grouped molecular tree and membership |
| `contrastive.nwk`, `.all.tsv`, `.sampled.tsv` | Raw selection, including unresolved candidates |
| `summary_tree.pdf`, `.svg` | Publication view: scientific names, group species counts, and red pair clades |
| `summary_tree_id.pdf`, `.svg` | Same tree with run, group, and pair IDs in tip labels |
| `phylogenetic_block_tree.pdf`, `.svg` | Same tree with expanded-clade shading and blue expanded branches; original pair branches remain red |
| `phylogenetic_block_tree_id.pdf`, `.svg` | Expanded view with run, group, original pair, and expanded-clade IDs |
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
From v0.8.13, IDs are ordered by clade membership rather than display sample
names. Recomputing an earlier result can therefore change its ID numbers;
reading or exporting saved results preserves their existing IDs.

## Phylogenetic blocks

`phylogenetic_block_id` assigns the largest rooted clade containing one selected
contrast pair and no other selected pair. Start at each pair's common ancestor
on the input molecular tree, then ascend until the parent contains another
pair. Unpaired tips remain available during expansion. The exported
`summary_tree.nwk` shows the same assignments using one display tip per trait
clade; the pair-only `contrastive.nwk` omits unpaired tips.

Blocks are disjoint and use the anchoring pair IDs. An added sample
can have a block ID and a blank original pair ID. Missing-trait
samples and branches outside the expanded clades have blank IDs. With a single
pair, its block can include the entire observed tree.

For leave-one-group-out validation, review the number of blocks and their
sample counts and label balance. A block can be much larger than
its anchor pair. See [PhenoRadar validation groups](phenoradar_inputs.md#validation-groups).

### Apply to saved results

Existing completed analyses are immutable. Export a derivative from saved pairs
without inference or pair reselection:

```bash
python workflow/scripts/expanded_pair_clades.py --analysis results/leaf/downstream/carnivory
```

Run with the workflow test environment or the pinned `timetree` environment.
Every configured contrast branch must already have its summary tree, metadata,
pair table, and summary. Outputs go under `<analysis>/phylogenetic_blocks/`:

- `metadata/species_metadata.tsv`: the full cohort, taxonomy, original pairs,
  and expanded-clade columns; its JSON records the joins.
- `phylogeny/<tree>/contrast/`: the unchanged saved tree and pair table,
  augmented branch metadata, expanded-clade counts, and both expanded figures.
- `provenance.json`: source/output checksums, code version, and branch counts.

Expression inputs remain at their original analysis paths. Publication is
atomic and existing export directories are refused; use `--outdir NEW_DIRECTORY`
for another export. Original files and completion records are unchanged.
