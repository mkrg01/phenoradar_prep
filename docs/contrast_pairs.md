# Trait contrast pairs

[Documentation](index.md) · [BUSCO phylogeny](phylogeny.md)

Pairs use nwkit's homogeneous-clade grouping and contrastive-clade selection.
`contrast.trait` (default `C4`) selects a column in `inputs.species_trait`;
missing traits stay unknown and expression replicates do not add species weight.

## Pairs from full or phenotyped trees

Choose `phylogeny.species_sets: [all]`, `[phenotyped]`, or both:

```bash
./run_pipeline.sh --cores 32 --resources mem_gb=128 -- phylogeny_contrast_pairs
```

This assigns pairs directly on the molecular tree and can schedule missing or
outdated inference. Pair-only jobs need 1 core/8 GB; check `--dry-run` before
reducing resources. `contrast.enabled` does not control this target.

Missing-trait tips are pruned with root direction/path lengths preserved in
`observed_tree.nwk`; their metadata retains empty pair IDs. Zero/one observed
state yields zero pairs; more than two states is unsupported. Outputs are in
`results/<run_name>/phylogeny/<set>/contrast/`.

## Representative analysis

```bash
./run_pipeline.sh --cores 32 --resources mem_gb=128 -- contrast_pairs
```

This needs [phylogeny inputs](phylogeny.md#inputs), exactly two observed states,
and at least four representatives. `contrast.enabled: true` includes it in
`all`; outputs go to `phylogeny/representatives/`.

1. Skim homogeneous clades on the local NCBI guide, selecting representatives
   by BUSCO completeness; top-level `seed` controls randomized ties.
2. Resolve an outgroup and infer a BUSCO tree for that subset.
3. Skim the molecular tree and pair minimal mixed clades with opposite-state
   representatives; multiway cases remain unresolved.
4. Map original species through both grouping stages to representatives/pairs.

A manual outgroup must belong to the representatives. Dating is not included.
Selection records are in `selection/`.

## After species exclusion

With nonempty `exclude_species`, `phylogeny_contrast_pairs` requires completed
[filtered results](species_filter.md) for the requested sets and never schedules
inference. Filtering recomputes pairs for completed full/phenotyped trees,
regardless of `phylogeny.species_sets`, under
`filtered/phylogeny/<set>/contrast/`. Representatives are not filtered.

## Outputs

Each `contrast/` directory contains:

| Output | Contents |
| --- | --- |
| `contrast_pairs.tsv` | Pair states, representatives, and species counts |
| `species_metadata.tsv` | Traits, groups, representatives, and nullable pair IDs |
| `summary_tree.nwk`, `.all.tsv`, `.sampled.tsv` | Molecular skim and membership |
| `contrastive.nwk`, `.all.tsv`, `.sampled.tsv` | Raw selection, including unresolved candidates |
| `summary_tree.pdf`, `summary_tree.svg` | Trait-colored tree with group/pair IDs |
| `summary.json` | Counts, settings, and provenance |

Zero-pair results are valid and produce header-only pair tables.

## Interpretation

Pairs describe sampled trait contrasts, not independent evolutionary origins.
NCBI-inherited membership is not molecularly tested for every species, and
assignments have no branch-support filter. Pair IDs belong to each result and
can change with species membership or exclusions.
