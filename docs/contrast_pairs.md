# Trait contrast pairs

[Documentation](index.md) · [BUSCO phylogeny](phylogeny.md)

Pairs use nwkit's homogeneous-clade grouping and contrastive-clade selection.
Top-level `trait` selects a column in `inputs.species_trait`; missing traits stay unknown.

## Configuration and execution

Set these options before [preparing analysis](datasets.md#run-an-analysis):

```yaml
trait: carnivory
phylogeny:
  trees: [representatives]
  contrast_pairs:
    enabled: true
```

```bash
./run_analysis.sh submit --analysis analyses/analysis001 --target contrast_pairs
```

The same target works for `trees: [all]`, `[phenotyped]`, `[representatives]`, or
several trees. It schedules missing prerequisites for exactly those species sets,
then selects pairs on each resulting tree. The enabled flag is required for the
explicit target and includes pairs in `all`. Use `phylogeny` to stop at inference.

| Tree | When species are reduced |
| --- | --- |
| `all` | Infer all selected species, then prune missing-trait tips for pairs |
| `phenotyped` | Remove missing-trait species before inference |
| `representatives` | Remove missing-trait species and choose NCBI representatives before inference |

For full/phenotyped trees, missing-trait tips are pruned with root direction/path
lengths preserved in `observed_tree.nwk`; their metadata retains empty pair IDs.
Zero/one observed state yields zero pairs; more than two states is unsupported.
Pair jobs default to 1 CPU/8 GB; missing inference steps use their own resources.

## Representative selection

Representative inference needs [phylogeny inputs](phylogeny.md#inputs), exactly
two observed states, and at least four representatives. It is available with
`trees: [representatives]` even when pair selection is disabled.

1. Skim homogeneous trait clades on the local NCBI guide, selecting representatives
   by BUSCO completeness; top-level `seed` controls randomized ties.
2. Resolve an outgroup and infer a BUSCO tree for that subset.
3. If pairs are enabled, skim the molecular tree and pair minimal mixed clades
   with opposite-state representatives; multiway cases remain unresolved.
4. Map original species through both grouping stages to representatives/pairs.

Both trait states are compressed. A manual outgroup must belong to the
representatives. Selection records are in `phylogeny/representatives/selection/`.
Dating and taxonomy checks currently do not support representative trees.

## After species exclusion

Normal `analysis.yaml` exclusions apply before inference for every selected tree,
including representatives. The separate [post hoc export](species_filter.md)
recomputes pairs only for completed `all`/`phenotyped` trees; it cannot export
representative results or replace inference on a newly selected species set.

## Outputs

Each `results/<analysis>/phylogeny/<tree>/contrast/` directory contains:

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
