# BUSCO phylogeny

[Documentation](index.md) · [Dating](dating.md) · [Contrast pairs](contrast_pairs.md)

Infer rooted trees from BUSCO markers, with one tip per analysis sample:

```text
BUSCO CDS → cdskit → FAMSA → trimAl → VeryFastTree → ASTRAL-IV/CASTLES-II
```

Branch lengths are substitutions/site; absolute ages require [dating](dating.md).
Multiple samples of one biological species remain separate tips and do not
represent additional evolutionary origins. Inference does not use OG/KO assignments.

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

For an input audit, marker plan, and outgroup check before inference, submit
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

`auto` chooses a single-species basal lineage using NCBI taxonomy, with nwkit's
APG IV order tree as an angiosperm fallback. For a basal group containing several
species, provide an explicit list. Review `rooting/outgroup.json`.

ASTRAL uses one sample as its root anchor; the output is then rooted on the
complete basal split. Taxonomy does not constrain the molecular topology.
The substitution tree splits the root edge at its midpoint; LSD2 re-estimates
that position. ASTRAL's single-anchor branch estimation for a multi-tip basal
clade has not been independently benchmarked here; provenance records this limitation.

## Marker and sequence selection

Only unambiguous single-copy `Complete` hits are used. Duplicated, fragmented,
multiply reported hits and genes assigned to several markers are omitted.
Coverage counts distinct biological species, so extra samples do not inflate it.

| Setting under `phylogeny` | Default | Meaning |
| --- | --- | --- |
| `max_markers` | `500` | Select by species coverage, then BUSCO ID |
| `min_taxa` | `4` | Minimum species per marker before and after QC; at least four |
| `min_protein_length` | `100` | Minimum known amino acids per extracted/trimmed sequence |
| `max_unknown_fraction` | `0.05` | Maximum unknown fraction in prepared proteins |
| `trimal_mode` | `gappyout` | `gappyout` or `automated1` |

cdskit pads, masks, and translates CDS using the build's genetic code. Stops and
unresolved codons become X. Review `species/*.json`: padding may affect reading
frames, and QC cannot verify ORFs or exclude paralogy.

FAMSA aligns proteins; trimAl selects columns with X treated as gaps. Retained
loci must pass sequence-length/species-coverage checks and contain a variable
amino-acid site. Markers lost during QC are not replaced.

## Inference and outputs

VeryFastTree uses double precision and `-lg -gamma` (LG+CAT search with Gamma20
length rescaling), reporting SH-like local support. ASTRAL-IV reports local
posterior probabilities and CASTLES-II substitution lengths. Both use seed
`12345`; no bootstrap or support filter is applied.

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

Marker alignment and gene-tree jobs default to 4 CPUs/8 GB; species-tree inference
uses 32 CPUs/64 GB. See [resource overrides](running.md#resource-budgets).
Logs are under `results/<build>/logs/downstream/<analysis>/phylogeny/`.

## Methods

[BUSCO](https://busco.ezlab.org/busco_userguide.html) ·
[cdskit](https://github.com/kfuku52/cdskit/tree/0.27.0/cdskit) ·
[FAMSA](https://github.com/refresh-bio/FAMSA) ·
[trimAl](https://github.com/inab/trimal) ·
[VeryFastTree](https://github.com/citiususc/veryfasttree) ·
[FastTree models/support](https://morgannprice.github.io/fasttree/) ·
[ASTRAL-IV/CASTLES-II](https://github.com/chaoszhang/ASTER/blob/master/tutorial/astral4.md)
