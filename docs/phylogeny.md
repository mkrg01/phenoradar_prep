# BUSCO species trees

[Documentation](index.md) · [Dating](dating.md)

The `phylogeny` target infers rooted species trees from existing BUSCO full tables
and original CDS, independently of ODB/KEGG:

```text
BUSCO markers -> cdskit -> FAMSA -> trimAl/QC -> VeryFastTree -> ASTRAL-IV/CASTLES-II
```

Branch lengths are substitutions/site; absolute ages require [dating](dating.md).

## Inputs

Alongside ordinary [inputs](inputs.md), provide one full table per species in
`phylogeny.busco_full_dir` (default `input/busco/full`). Use one lineage
dataset/version; `phylogeny.lineage` defaults to `embryophyta_odb12`.

Full tables need `Busco id`, `Status`, `Sequence`, `Score`, and `Length` columns,
a lineage header, and the same complete marker-ID set. Counts-only summaries and
genome-mode tables are unsupported. Supported paths below also accept `.gz`;
missing or multiple matches are errors:

```text
{species}.busco.full.tsv
{species}.tsv
{species}/full_table.tsv
{species}/run_{lineage}/full_table.tsv
```

Supply original, oriented, in-frame CDS matching BUSCO hits through `inputs.cds_dir`.
IDs such as `Species_g123:60-698` resolve to `Species_g123`; the full CDS is
translated, without reconstructing BUSCO-predicted peptides. Missing or ambiguous
original IDs fail extraction.

## Species sets

```yaml
phylogeny:
  species_sets: [all, phenotyped]
  trait: C4
```

| Set | Species | Output under `results/<run_name>/` |
| --- | --- | --- |
| `all` (default) | All species passing input selection | `phylogeny/all/` |
| `phenotyped` | Selected species with a nonmissing `phylogeny.trait` | `phylogeny/phenotyped/` |

Each set needs `min_taxa` species and gets independent markers, alignments,
trees, and roots. The `all` set needs no trait file; see [traits](inputs.md#traits)
for defining the phenotyped subset.

Dating, taxonomy checks, and `phylogeny_contrast_pairs` also follow `species_sets`.
[Representative analysis](contrast_pairs.md#representative-analysis) is separate.

## Setup and execution

```bash
# Optional input audit, marker plan, and outgroup check.
./run_pipeline.sh --cores 4 --resources mem_gb=16 -- phylogeny_prepare

# Infer gene trees and species trees.
./run_pipeline.sh --cores 32 --resources mem_gb=128 -- phylogeny
```

Set `phylogeny.enabled: true` to include trees in `all`. Tools are bundled in the
container; see [native execution](containers.md#native-execution) for Conda setup.

## Rooting

`phylogeny.outgroup` accepts `auto` or an exact species ID in every requested set.
Automatic selection uses NCBI taxonomy, with nwkit's APG IV order tree as an
angiosperm fallback. It chooses within each set; supply an outgroup if it fails.

Review `rooting/outgroup.json`. The outgroup is supplied before CASTLES-II
length estimation; the taxonomy guide does not constrain molecular topology.

## Marker and sequence selection

Only unambiguous single-copy `Complete` hits are eligible. Duplicated, fragmented,
missing, multiply reported hits, and genes assigned to multiple markers are
omitted. Coverage counts species, not expression runs.

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
species and a variable amino-acid site. Changing trimAl mode reuses raw alignments.

VeryFastTree uses double precision and `-lg -gamma` (LG+CAT search, Gamma20 length
rescaling), with SH-like local supports and no bootstrap or support filtering.
ASTRAL-IV combines gene trees, reporting local posterior probabilities and
CASTLES-II substitution lengths. Both use the top-level `seed`.

Every species must occur in a retained gene tree; review `species_coverage.tsv`.
Missing data, gene-tree error, paralogy, and model assumptions affect estimates.

## Resources

Per-job memory totals all threads; see [overrides](running.md#resource-budgets).

| Rule | Job unit | Threads | Memory (GB) |
| --- | --- | --- | --- |
| `align_busco_marker`, `infer_busco_gene_tree` | Marker | 4 | 8 |
| `infer_busco_species_tree` | Species set | 32 | 64 |

Preparation, extraction, collection, trimming, merging, calibration retrieval,
and dating each default to 1 thread and 4 GB per job.

## Outputs

Under `results/<run_name>/phylogeny/<set>/`:

| Output | Contents |
| --- | --- |
| `selection/`, `rooting/` | Phenotyped selection, outgroup, and supporting records |
| `plan/` | Marker ranks/selection, source paths, lineage, and provenance |
| `species/*.faa`, `*.json` | Prepared proteins and sequence QC |
| `markers/`, `alignments/raw/` | Marker inputs and raw alignments |
| `alignments/*.faa`, `*.json`, `*.columns.tsv` | Trimmed alignments, QC, and 1-based final-to-raw columns |
| `gene_trees/`, `gene_trees.nwk`, `gene_trees.json` | Individual/merged gene trees and diagnostics |
| `species_coverage.tsv` | Retained locus counts per species |
| `species_tree.nwk`, `species_tree.json` | Rooted species tree and inference provenance |

Logs/benchmarks are under `logs/<run_name>/phylogeny/`. Optional dating, taxonomy
review, and contrast-pair outputs live beside their source tree.

## Methods and source documentation

[BUSCO formats](https://busco.ezlab.org/busco_userguide.html) ·
[cdskit](https://github.com/kfuku52/cdskit/tree/0.27.0/cdskit) ·
[FAMSA](https://github.com/refresh-bio/FAMSA) ·
[trimAl](https://github.com/inab/trimal) ·
[VeryFastTree](https://github.com/citiususc/veryfasttree) ·
[FastTree models/support](https://morgannprice.github.io/fasttree/) ·
[ASTRAL-IV/CASTLES-II](https://github.com/chaoszhang/ASTER/blob/master/tutorial/astral4.md)
