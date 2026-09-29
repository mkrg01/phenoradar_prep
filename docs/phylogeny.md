# BUSCO trees of separately assembled samples

[Documentation](index.md) · [Dating](dating.md)

`phylogeny` infers rooted trees from the completed build's BUSCO full tables and
CDS, with one tip per analysis sample. Multiple samples of one biological species
remain separate tips; they are not additional species or independent evolutionary
origins. Tree inference does not use ODB/KEGG assignments:

```text
BUSCO markers -> cdskit -> FAMSA -> trimAl/QC -> VeryFastTree -> ASTRAL-IV/CASTLES-II
```

Branch lengths are substitutions/site; absolute ages require [dating](dating.md).

## Inputs

CDS and full tables are staged automatically from `results/<build>/database/`.
Analysis inherits `busco.lineage` from build settings (default `embryophyta_odb12`).

Full tables need `Busco id`, `Status`, `Sequence`, `Score`, and `Length` columns,
a lineage header, and the same complete marker-ID set. Counts-only summaries and
genome-mode tables are unsupported. Complete the build before starting analysis.

Original CDS must be oriented, in-frame, and match the BUSCO hits.
IDs such as `Species_accession_g123:60-698` resolve to `Species_accession_g123`; the full CDS is
translated, without reconstructing BUSCO-predicted peptides. Missing or ambiguous
original IDs fail extraction.

## Species sets

Set options before [preparing analysis](datasets.md#run-an-analysis):

```yaml
trait: C4
phylogeny:
  trees: [representatives]
  contrast_pairs:
    enabled: false
```

| Tree | Inference samples | Output under `results/<build>/downstream/<analysis>/` |
| --- | --- | --- |
| `all` | All samples passing input selection | `phylogeny/all/` |
| `phenotyped` | Selected samples whose biological species has a nonmissing `trait` | `phylogeny/phenotyped/` |
| `representatives` | Known-trait samples compressed on the NCBI guide | `phylogeny/representatives/` |

Use `trees: []` to disable inference, or list multiple trees to infer them
independently. Each tree gets its own markers, alignments, roots, and diagnostics.
Representative selection requires two observed states and at least four
representatives; both trait states are compressed, selecting by BUSCO completeness
with seeded ties. See [representative selection](contrast_pairs.md#representative-selection).
The `all` tree needs no trait file. See [traits](inputs.md#traits) for missing values.

Pair selection, dating, and taxonomy checks use the selected trees and require
a nonempty `trees` list. Dating and taxonomy checks reject representative trees.

## Setup and execution

```bash
# Optional input audit, marker plan, and outgroup check.
./run_analysis.sh submit --analysis results/angiosperm_leaf_20260925/downstream/c4_photosynthesis_20260929 --target phylogeny_prepare

# After that job finishes and the audit is reviewed:
./run_analysis.sh submit --analysis results/angiosperm_leaf_20260925/downstream/c4_photosynthesis_20260929 --target phylogeny
```

`phylogeny` includes missing preparation steps and stops at inference.
`all` includes selected trees and enabled postprocessing.

## Rooting

`phylogeny.outgroup` accepts `auto`, a biological species ID, a sample ID, or a
list such as `[Basal_species_A, Basal_species_B]`. An explicit species/sample
selector expands to every selected sample of that biological species. The whole
basal group must form one side of a split in the molecular tree. All samples are
retained; a nonmonophyletic group stops inference instead of pruning taxa.

Automatic selection conservatively chooses a complete single-species basal
lineage using NCBI taxonomy, with nwkit's APG IV order tree as an angiosperm
fallback. It retains all samples of that species. Specify a list explicitly for
a basal group containing several species. No extra species are added.

Review `rooting/outgroup.json`. ASTRAL-IV accepts only one root label: the first
sample ID in the sorted basal group is used as a reproducible anchor during
CASTLES-II estimation. The resulting tree is then oriented on the complete
basal split, preserving tips and unrooted edges. Numeric supports follow their
edges. The taxonomy guide does not constrain molecular topology.

Root placement along the selected edge is a midpoint convention in the
substitution tree. [LSD2 dating](dating.md) re-estimates that placement while
keeping the edge and topology fixed. The effect of ASTRAL's single-anchor
branch estimation with a multi-tip basal clade has not been independently
benchmarked here; this limitation is recorded in species-tree provenance.

## Marker and sequence selection

Only unambiguous single-copy `Complete` hits are eligible. Duplicated, fragmented,
missing, multiply reported hits, and genes assigned to multiple markers are
omitted. Each eligible sample contributes its own marker sequence and tree tip. Coverage
and `min_taxa` count distinct biological `species_id` values, including after
alignment QC; extra samples do not inflate species coverage. NCBI guide trees
expand a shared taxid to multiple sample tips without assigning biological branch lengths.

| Setting under `phylogeny` | Default | Use |
| --- | --- | --- |
| `max_markers` | `500` | Highest coverage first, then BUSCO ID ascending |
| `min_taxa` | `4` | Minimum species per marker before and after QC; at least four |
| `min_protein_length` | `100` | Known amino acids per extracted/trimmed sequence |
| `max_unknown_fraction` | `0.05` | Maximum unknown fraction in prepared proteins |
| `trimal_mode` | `gappyout` | `gappyout` or `automated1` |

Markers lost during QC are not replaced. cdskit `pad`, `mask`, and `translate`
use `translation.table`; stops/unresolved codons become X. Review `species/*.json`:
padding can change reading frames, and QC cannot verify ORFs or exclude paralogy.

## Alignment and tree inference

FAMSA aligns proteins; trimAl selects columns with X treated as gaps, then restores
X in retained columns. After trimming, sequences need `min_protein_length` known
residues; columns without known residues are removed. Each locus needs `min_taxa`
species and a variable amino-acid site.

VeryFastTree uses double precision and `-lg -gamma` (LG+CAT search, Gamma20 length
rescaling), with SH-like local supports and no bootstrap or support filtering.
ASTRAL-IV combines gene trees, reporting local posterior probabilities and
CASTLES-II substitution lengths. Both use the top-level `seed`.

Every selected sample must occur in a retained gene tree; review `species_coverage.tsv`.
Missing data, gene-tree error, paralogy, and model assumptions affect estimates.

## Resources

Per-job defaults; see [resource overrides](running.md#resource-budgets).

| Rule | Job unit | Threads | Memory (GB) |
| --- | --- | --- | --- |
| `align_busco_marker`, `infer_busco_gene_tree` | Marker | 4 | 8 |
| `infer_busco_species_tree` | Species set | 32 | 64 |

Preparation/QC steps generally use 1 CPU/4 GB; [dating](dating.md#resources) requires 1 CPU.

## Outputs

Under `results/<build>/downstream/<analysis>/phylogeny/<set>/`:

| Output | Contents |
| --- | --- |
| `selection/`, `rooting/` | Phenotyped/representative selection, outgroup, and supporting records |
| `plan/` | Marker ranks/selection, source paths, lineage, and provenance |
| `species/*.faa`, `*.json` | Prepared proteins and sequence QC |
| `markers/`, `alignments/raw/` | Marker inputs and raw alignments |
| `alignments/*.faa`, `*.json`, `*.columns.tsv` | Trimmed alignments, QC, and 1-based final-to-raw columns |
| `gene_trees/`, `gene_trees.nwk`, `gene_trees.json` | Individual/merged gene trees and diagnostics |
| `species_coverage.tsv` | Retained locus counts per analysis sample |
| `species_tree.nwk`, `species_tree.json` | Rooted species tree and inference provenance |

Logs/benchmarks are under `results/<build>/logs/downstream/<analysis>/phylogeny/`. Optional dating, taxonomy
review, and contrast-pair outputs live beside their source tree.

## Methods and source documentation

[BUSCO formats](https://busco.ezlab.org/busco_userguide.html) ·
[cdskit](https://github.com/kfuku52/cdskit/tree/0.27.0/cdskit) ·
[FAMSA](https://github.com/refresh-bio/FAMSA) ·
[trimAl](https://github.com/inab/trimal) ·
[VeryFastTree](https://github.com/citiususc/veryfasttree) ·
[FastTree models/support](https://morgannprice.github.io/fasttree/) ·
[ASTRAL-IV/CASTLES-II](https://github.com/chaoszhang/ASTER/blob/master/tutorial/astral4.md)
