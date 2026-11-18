# BUSCO phylogeny

[Documentation](index.md) · [Dating](dating.md) · [Contrast pairs](contrast_pairs.md)

Infer rooted trees from BUSCO markers, with one tip per analysis sample:

```text
BUSCO CDS → cdskit → FAMSA → trimAl → VeryFastTree → ASTRAL-IV/CASTLES-II
```

Branch lengths are in substitutions per site. For ages, see [dating](dating.md).
Multiple samples of a biological species remain separate tips.

## Species sets

Choose trees in `analysis.yaml` before [submitting an analysis](datasets.md#run-an-analysis):

```yaml
phylogeny:
  trees: [all]
  outgroup: auto
```

| Tree | Samples used |
| --- | --- |
| `all` | All samples passing analysis selection; no traits required |
| `phenotyped` | Selected samples with a known species trait |
| `representatives` | Known-trait samples reduced using the NCBI guide and BUSCO completeness |

List several trees to infer them independently, or `[]` to disable inference.
Representative selection needs two observed trait states and at least four
representatives; see [selection details](contrast_pairs.md#representative-selection).
Dating and taxonomy checks support only `all` and `phenotyped`.

## Setup and execution

```bash
./run_analysis.sh submit --analysis results/leaf/downstream/carnivory --target phylogeny
```

For an input audit, marker plan, and ASTRAL anchor check before inference, submit
`--target phylogeny_prepare` first. `phylogeny` stops at inference; `all` also
runs enabled postprocessing.

## Inputs

CDS and BUSCO full tables are staged from the completed database. Full tables
must contain `Busco id`, `Status`, `Sequence`, `Score`, and `Length`, a lineage
header, and the same complete marker-ID set. Counts-only summaries and genome-mode
tables are unsupported.

Original CDS must be oriented, in-frame, and match BUSCO hits. A hit such as
`Species_run_g123:60-698` resolves to `Species_run_g123`; the full original CDS
is translated. Missing or ambiguous IDs stop extraction.

## Rooting

`phylogeny.outgroup` accepts `auto`, a species/sample ID, or a list of IDs.
Each selector expands to all selected samples of that biological species.
The complete group must lie on one side of a split in the molecular tree;
a nonmonophyletic group stops inference. No extra species are added.

`auto` chooses the smallest basal lineage in the frozen NCBI guide by species
count and uses its most BUSCO-complete sample as the ASTRAL anchor. An explicit
outgroup uses the first sample in sorted order. See `rooting/anchor.json`.

After inference, nwkit tries NCBI, OpenTree, and TimeTree in that order. The NCBI
step uses the frozen taxonomy. The other two query their services by scientific
name. If the references cannot resolve a split, check `rooting/outgroup.log`
under the analysis logs and specify `phylogeny.outgroup` explicitly. Review
`rooting/outgroup.json` for the final split.

Rooting preserves the inferred topology and branch lengths. A new root splits
its edge at the midpoint. treePL uses the rooted split for dating. The branch lengths
retain the effect of ASTRAL's single internal anchor. The ASTRAL output and
provenance before final rooting are in `rooting/astral.nwk` and
`rooting/astral.json`.

## Marker and sequence selection

Only unambiguous single-copy `Complete` hits are used. Duplicated, fragmented,
multiply reported hits and genes assigned to several markers are omitted.
Coverage counts distinct biological species, so extra samples do not inflate it.

| Setting under `phylogeny` | Default | Meaning |
| --- | --- | --- |
| `max_markers` | `null` | All eligible markers; a positive integer caps by species coverage, then BUSCO ID |
| `min_taxa` | `0.5` | Minimum fraction (0 to 1) of the tree's selected species represented per marker |
| `min_protein_length` | `100` | Minimum known amino acids per extracted/trimmed sequence |
| `max_unknown_fraction` | `0.05` | Maximum unknown fraction in prepared proteins |
| `trimal_mode` | `gappyout` | `gappyout` or `automated1` |

cdskit pads, masks, and translates CDS using the build's genetic code. Stops and
unresolved codons become X. Review `species/*.json`: padding may affect reading
frames, and QC cannot verify ORFs or exclude paralogy.

FAMSA aligns proteins and trimAl selects columns with X treated as gaps. The
columns selected by trimAl are preserved when short sequences are removed;
columns left entirely as gaps or X are not removed afterward. The
same species coverage threshold applies at marker selection, before alignment,
and after trimming and sequence filtering. Each check uses all species selected
for that tree as the denominator. No separate four-species or variable-site
filter is applied to markers. Variable and parsimony-informative site counts
are recorded for diagnosis. Gene trees require at least two tips for CASTLES-II
branch-length estimation; the complete analysis still requires at least four
species. The marker list is fixed at selection.

## BUSCO alignments

Raw and trimmed marker alignments are in `phylogeny/<set>/alignments/`.
BUSCO uses single-copy hits from the samples in the chosen tree.
[OG alignments](alignments.md) are a separate output based on OrthoDB and
include all mapped copies from the analysis cohort.

## Inference and outputs

VeryFastTree uses double precision and `-lg -gamma` (LG+CAT search with Gamma20
length rescaling), reporting SH-like local support. ASTRAL-IV reports local
posterior probabilities and CASTLES-II substitution lengths. Both use seed
`12345`; no bootstrap or support filter is applied.
VeryFastTree defaults to one thread: repeated four-thread runs of the pinned
4.0.5 binary produced different topologies despite the same input and seed.

Under `results/<build>/downstream/<analysis>/phylogeny/<set>/`:

| Output | Review |
| --- | --- |
| `species_tree.nwk`, `species_tree.json` | Rooted tree and inference provenance |
| `species_coverage.tsv` | Retained loci per sample; every selected sample must appear |
| `selection/`, `rooting/`, `plan/` | Sample selection, outgroup, and marker choices |
| `species/*.json`, `alignments/*.json` | Sequence and alignment QC |
| `alignments/*.faa`, `*.columns.tsv` | Trimmed alignments and final-to-raw column positions |
| `gene_trees/`, `gene_trees.nwk`, `gene_trees.json` | Locus trees and diagnostics |

Review coverage and rooting before interpreting the tree. Missing data,
gene-tree error, paralogy, and model assumptions affect the estimates.

## Resources

On Slurm, each BUSCO marker's alignment, trimming, and gene-tree inference
share one worker job. FAMSA alignment requests 1 CPU and 16 GB for each tree
set, so a full group defaults to 1 CPU and 16 GB. Species-tree inference is
a separate job with a default of 32 CPUs and 128 GB. See
[resource overrides](running.md#resource-budgets). Logs are under
`results/<build>/logs/downstream/<analysis>/phylogeny/`.

## Methods

[BUSCO](https://busco.ezlab.org/busco_userguide.html) ·
[cdskit](https://github.com/kfuku52/cdskit/tree/0.27.0/cdskit) ·
[FAMSA](https://github.com/refresh-bio/FAMSA) ·
[trimAl](https://github.com/inab/trimal) ·
[VeryFastTree](https://github.com/citiususc/veryfasttree) ·
[FastTree models/support](https://morgannprice.github.io/fasttree/) ·
[ASTRAL-IV/CASTLES-II](https://github.com/chaoszhang/ASTER/blob/master/tutorial/astral4.md)
