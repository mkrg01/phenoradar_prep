# Contrast pairs from a representative species tree

[Back to README](../README.md)

The optional `contrast_pairs` target runs independently of the full species
tree, ODB/KEGG mapping and expression aggregation. It uses the same BUSCO
inference rules, tool environments, marker cap (500 by default) and sequence
QC as `phylogeny`, with a separate representative manifest and output directory.
It does not run dating. Full and subset alignments/gene trees are separate
because marker coverage and alignments depend on the selected species set.

## Inputs and execution

Use the ordinary metadata/BUSCO/CDS inputs and `phylogeny.busco_full_dir`.
Species selection currently also validates the usual abundance input files,
although contrast inference does not compute expression tables.
All phenotype annotations come from `inputs.species_trait`, never metadata:

```text
species	C4
Plant alpha	0
Plant beta	1
Plant gamma	
```

Space-separated names and normalized underscore IDs are accepted. Hyphens
are preserved, as in the existing phylogeny labels. Duplicate normalized names
and invalid C4 values fail validation. Blank/NA annotations and species absent
from the trait table remain unknown. Extra species in that table do not enter
the analysis. The observed species must have exactly two states.

```yaml
inputs:
  species_trait: species_trait/species_trait.tsv
phylogeny:
  busco_full_dir: /path/to/busco_full_longest_cds
  sequence_dir: /path/to/longest_cds
  outgroup: auto
contrast:
  enabled: false
  trait: C4
```

`phylogeny.outgroup: auto` enables
[automatic outgroup selection](phylogeny.md#setup-and-execution). First skim
selects the inference species; the outgroup is then chosen **within that exact
representative set**. No additional species is included for rooting. The full
and contrast branches apply the same method to their own species sets and can
choose different outgroups. A manual species label must already occur in the
relevant inference manifest, including the skim representatives for contrast.
The chosen outgroup retains its trait/group information and participates in
the subsequent skims. The figure uses ordinary species labels; the rooting
record identifies which representative was used as the outgroup.

From the project root, after building ASTRAL as described in the phylogeny guide:

```bash
./run_pipeline.sh --software-deployment-method conda \
  --configfile config/mydata.yaml config/phylogeny.local.yaml \
  --cores 16 --resources mem_gb=80 --dry-run -- contrast_pairs

./run_pipeline.sh --software-deployment-method conda \
  --configfile config/mydata.yaml config/phylogeny.local.yaml \
  --cores 16 --resources mem_gb=80 -- contrast_pairs
```

Set `contrast.enabled: true` to include these outputs in `all`. The explicit
target works while false. For preparation only, target
`results/<analysis>/contrast/selection/selection.json`. To also resolve the
outgroup without inferring trees, target
`results/<analysis>/contrast/phylogeny/rooting/outgroup.json`.

## Processing

1. Build the NCBI guide from the selected dataset.
   The pinned nwkit 0.27.0 lineage operations use the frozen SQLite snapshot;
   they need no traversal pickle or extra taxonomy download.
2. Restrict the guide to observed traits and run the first nwkit skim:
   one highest-completeness species per homogeneous clade. Completeness is
   `(single + duplicated) / total`, without display rounding. The input rows
   are sorted and `phylogeny.seed` fixes randomized ties.
3. Resolve the outgroup from the compressed NCBI guide, whose tips are exactly
   the representatives, and infer their BUSCO tree. At least four
   representatives are required. Marker coverage is
   ranked within this subset; the existing sequence/alignment QC remains active.
4. Skim the rooted molecular tree again, retaining the outgroup as a representative.
5. Apply nwkit's `only_contrastive_clades` selection. Minimal mixed clades with
   exactly two opposite-state representatives become pairs. Multiway cases
   remain unresolved; they are not arbitrarily decomposed into pairs.
6. Compose the first and second group maps to assign original species to final
   representatives and pair IDs. Export one summary figure in PDF and SVG.

The adapter calls the grouping, sampling and contrastive-clade functions used
by `nwkit skim`, without modifying nwkit. Raw nwkit group/contrastive IDs remain
in the per-stage tables. Public pair IDs are ordered by their representative
species names; IDs can change when membership changes between analyses.

## Outputs

Under `results/<analysis>/contrast/`:

| Output | Contents |
| --- | --- |
| `selection/ncbi_skim.nwk`, `.all.tsv`, `.sampled.tsv` | Initial compressed NCBI tree and complete first-stage membership |
| `selection/samples.tsv`, `traits.tsv`, `selection.json` | Inference manifest, dataset-wide normalized traits/roles, source hash and selection counts |
| `phylogeny/` | The same marker plans, alignments, gene trees, coverage and rooted species tree as the full branch |
| `phylogeny/rooting/outgroup.txt`, `outgroup.json` | Outgroup selected within the representative set, candidate manifest/count and selection evidence |
| `summary_tree.nwk`, `.all.tsv`, `.sampled.tsv` | Second skim on the inferred molecular tree |
| `contrastive.nwk`, `.all.tsv`, `.sampled.tsv` | Raw nwkit contrastive-clade selection, including unresolved multiway candidates |
| `contrast_pairs.tsv` | One row per pair: state values, representatives, final groups and original species counts |
| `species_metadata.tsv` | All selected species: original trait, role, final group/representative and nullable pair ID |
| `summary.json` | Pair counts, unresolved clades, source records and assignment interpretation |
| `summary_tree.pdf`, `summary_tree.svg` | Trait-colored tips, italic species names, and aligned columns for species counts (`n`), group IDs and bold pair IDs |

The figure uses Matplotlib already supplied with nwkit; no R environment is
required. Species names use 8 pt italic type; annotation columns use upright
type. The nominal width is 7.2 inches (183 mm), expanding only for unusually
long labels, with height adjusted to the tip count. The PDF embeds TrueType
fonts and SVG retains editable text. A scale bar reports substitutions/site;
the figure has no pair-count title or outgroup annotation. Branch lengths are
not dates. Ladderization only changes the order of siblings in the drawing.
Original ASTRAL
supports remain in the inferred and compressed molecular trees; support-based
pair filtering is not applied. Figure-only recovery reuses the inferred tree.

Unknown traits have no group/pair assignment. The outgroup is an observed
representative and is processed by the same grouping/pair criteria as other tips.
A species can be a pair member without being a representative. Pair membership
for non-representatives is **inherited from the NCBI group**, not separately
validated by molecular inference. These pairs are comparisons within the sampled
known-trait tree, not estimates of independent C4 origins. A zero-pair analysis
is successful with a header-only pair table. When nwkit finds no contrastive
tips at all, `contrastive.nwk` is empty and `summary.json` reports the result.

## Validation

`tests/test_contrast_pairs.py` covers the frozen NCBI guide, automatic root
selection, rejection of incomplete/multi-species root assignments, trait-name
normalization, missing traits, deterministic ties, composition of both group
maps, multiway candidates, empty contrastive outputs, and PDF/SVG generation.
The real-tool workflow test runs both inference branches with different
automatically selected roots. It checks that the full-tree outgroup is not
added to the contrast representatives, unchanged reruns, recovery of a deleted
figure, and isolation of the full tree from a trait-only change. It uses
synthetic sequences and a local taxonomy fixture, not external services.

On the current dataset, the first skim selects 199 representatives from 2,047
known-trait species. Root selection chooses `Nymphaea_colorata` within those
199 representatives, while the full 5,586-species manifest selects
`Amborella_trichopoda`. The full 199-species contrast inference has not been run.
A small real-data check selects eight known-trait representatives from nine
input species, roots on `Nymphaea_colorata`, and completes inference with 20
markers, three contrast pairs, and both figure formats. This checks execution
and membership; it does not establish the accuracy of the inferred phylogeny.
