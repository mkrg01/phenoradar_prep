# Trait contrast pairs

[Documentation](index.md) · [BUSCO phylogeny](phylogeny.md)

Two targets assign pairs using nwkit's homogeneous-clade grouping and
contrastive-clade selection:

| Target | Tree used | Membership of non-representative species |
| --- | --- | --- |
| `phylogeny_contrast_pairs` | Full or phenotyped BUSCO tree | Determined directly on the molecular tree |
| `contrast_pairs` | Newly inferred tree of trait-guided NCBI representatives | Inherited through the NCBI and molecular group maps |

Both use `inputs.species_trait`, with `contrast.trait: C4` by default. See
[trait formats](inputs.md#traits). Missing traits remain unknown. Multiple RNA-seq
runs do not increase a species' weight. Pair IDs belong to one result directory
and can change when membership changes.

## Pairs from full or phenotyped trees

Set `phylogeny.species_sets` to `[all]`, `[phenotyped]`, or both, then run:

```bash
./run_pipeline.sh --configfile config/mydata.yaml \
  --cores 32 --resources mem_gb=128 \
  -- phylogeny_contrast_pairs
```

With `exclude_species: []`, this follows normal inference dependencies:
completed unchanged trees are reused, and missing or outdated trees can schedule
inference. Use the inference budget above unless a dry-run confirms that only
pair assignment remains. A budget of one core and 8 GB fits the declared
pair-only jobs.
`contrast.enabled` does not control this target.

| Input under `results/<run_name>/` | Pair results |
| --- | --- |
| `phylogeny/all/species_tree.nwk` | `phylogeny/all/contrast/` |
| `phylogeny/phenotyped/species_tree.nwk` | `phylogeny/phenotyped/contrast/` |

Inputs are the rooted tree and QC record, its sample manifest,
`metadata/metadata_high_busco.tsv` for representative scores, and the trait table.
`contrast.trait` may differ from the trait used to select phenotyped inference.

The tree must match its manifest and recorded outgroup. Missing-trait tips are
pruned with surviving path lengths and root direction preserved. The observed
subtree may therefore lack the original outgroup. `observed_tree.nwk` records
this subtree before grouping; no new outgroup is selected.

The molecular skim selects representatives and assigns pairs directly from this
tree. The output species metadata retains every species in the source inference
set, including unknown-trait species with empty group/pair IDs. Zero or one
observed state produces zero pairs; more than two states is unsupported. There
is no four-representative minimum for this postprocessing step. With no observed
species, Newick files are empty, TSVs retain headers, and figures explain the
empty result.

## Representative analysis

`contrast_pairs` first selects representatives on an NCBI guide, then infers their
BUSCO tree using the same sequence QC, marker cap, and tools as `phylogeny`.
It needs the ordinary metadata/BUSCO/CDS inputs, valid abundance paths during
selection, and `phylogeny.busco_full_dir`. The observed species must have exactly
two states and yield at least four representatives.

Run:

```bash
./run_pipeline.sh --configfile config/mydata.yaml \
  --cores 32 --resources mem_gb=128 -- contrast_pairs
```

Set `contrast.enabled: true` to include this analysis in `all`. Its output is
`results/<run_name>/phylogeny/representatives/`. Marker coverage is ranked within
that subset; alignments and gene trees are inferred independently. Dating is
not part of this target.

For selection alone, target
`results/<run_name>/phylogeny/representatives/selection/selection.json`.
To also resolve the outgroup, target
`results/<run_name>/phylogeny/representatives/rooting/outgroup.json`.
A manual `phylogeny.outgroup` must be among the selected representatives;
`auto` chooses within that set.

1. Build the NCBI guide from the selected dataset and frozen taxonomy snapshot.
2. Keep observed-trait species and skim homogeneous clades, choosing the species
   with the highest complete BUSCO fraction in each. Sorted inputs and
   `phylogeny.seed` make randomized ties reproducible.
3. Resolve the root within the representative set and infer its BUSCO tree.
4. Skim the rooted molecular tree, then apply nwkit's `only_contrastive_clades`.
   Minimal mixed clades with two opposite-state representatives form pairs;
   multiway cases remain unresolved.
5. Compose the NCBI and molecular group maps to assign original species to final
   representatives and pairs.

The outgroup participates in grouping and pair selection with its observed trait.
The adapter uses the pinned nwkit grouping and sampling functions; raw stage IDs
are retained alongside public pair IDs ordered by representative species names.

## After species exclusion

With a nonempty `exclude_species`, `phylogeny_contrast_pairs` uses the same
completed-result export as [filter_species](species_filter.md). It requires ready
inputs for the requested species sets and never schedules inference in this mode.
`filter_species` instead exports what is complete and reports unavailable sections
in `filtered/manifest.json`.

Filtering recomputes pairs for both completed molecular trees, independent of
`phylogeny.species_sets`, even if previous pair outputs do not exist. It needs
the tree/QC, manifest, BUSCO scores, and traits; raw sequences and gene trees are
unnecessary for pair assignment. Results go to `filtered/phylogeny/all/contrast/`
and `filtered/phylogeny/phenotyped/contrast/` when ready.

Exports use the original trees and root direction, even if the outgroup is excluded.
Changing exclusions can change pair IDs. The representative analysis is not filtered.

## Outputs

Both routes write the following in their branch's `contrast/` directory:

| File | Contents |
| --- | --- |
| `summary_tree.nwk`, `.all.tsv`, `.sampled.tsv` | Molecular skim tree and membership |
| `contrastive.nwk`, `.all.tsv`, `.sampled.tsv` | Raw contrastive-clade selection, including unresolved multiway candidates |
| `contrast_pairs.tsv` | Pair states, representatives, final groups, and original species counts |
| `species_metadata.tsv` | Source species, trait, role, final group/representative, and nullable pair ID |
| `summary.json` | Counts, unresolved clades, source checksums, settings, and assignment interpretation |
| `summary_tree.pdf`, `summary_tree.svg` | Trait-colored tree with species counts, group IDs, and pair IDs |

Full/phenotyped pair outputs also include `observed_tree.nwk`.
The representative route adds `selection/ncbi_skim.nwk`, `.all.tsv`, `.sampled.tsv`,
`samples.tsv`, `traits.tsv`, and `selection.json`, plus the normal BUSCO inference
and rooting outputs beside `contrast/`.

Figures retain editable text in SVG and embedded TrueType fonts in PDF. Their
scale bar uses substitutions/site. Ladderization changes drawing order only.
Original ASTRAL supports remain in molecular Newick files; pair assignment does
not filter by support.

A zero-pair result is successful and has a header-only pair table. When no
contrastive tips exist, `contrastive.nwk` is empty. Figure-only recovery reuses
pair assignment and inference. To recover a deleted file inside the filtered
bundle, use `--forcerun filter_species` with the `filter_species` target.

## Interpretation

Pairs describe comparisons in the sampled known-trait tree, not independent
origins of C4. A species can belong to a pair without being a representative.
In the representative route, membership inherited from an NCBI group has not
been separately tested by molecular inference for every member. Missing-trait
species have no assignment, and multiway contrasts are left unresolved.
