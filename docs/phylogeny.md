# BUSCO species trees

[Documentation](index.md) · [Dating](dating.md)

The `phylogeny` target infers a rooted species tree from existing BUSCO full
tables and original CDS or proteins. It uses a separate sequence preparation
path from ODB/KEGG and requires no orthogroup mapping or expression aggregation.

```text
BUSCO markers -> cdskit CDS preparation -> FAMSA -> trimAl -> alignment QC
  -> VeryFastTree gene trees -> ASTRAL-IV with CASTLES-II -> species_tree.nwk
```

The species tree has branch lengths in substitutions per site. Absolute ages
require the optional [dating stage](dating.md).

## Inputs

Use the ordinary [sample metadata and BUSCO summary](inputs.md), plus full
BUSCO tables in `phylogeny.busco_full_dir` (default `input/busco/full`).
Selected samples still need valid CDS and abundance paths during metadata
preparation, even when only tree inference is requested.

Full tables must have `Busco id`, `Status`, `Sequence`, `Score`, and `Length`
columns, a lineage header, and the same complete marker-ID set. Use one lineage
dataset/version throughout. A counts-only summary is insufficient; genome-mode
tables are unsupported because they do not provide enough information to
reconstruct spliced proteins.

| Setting under `phylogeny` | Default | Meaning |
| --- | --- | --- |
| `busco_full_dir` | `input/busco/full` | Full-table directory |
| `busco_full_suffix` | `.busco.full.tsv` | Suffix after exact species ID; a configured `.tsv.gz` suffix is supported |
| `lineage` | `embryophyta_odb12` | Expected BUSCO lineage |
| `sequence_mode` | `cds` | `cds` or `protein` |
| `sequence_dir` | `null` | Reuse CDS paths in `samples.tsv`; required in protein mode |
| `sequence_suffix` | `_longestCDS.fa.gz` | Suffix when using `sequence_dir` |

CDS inputs must be the original, oriented, in-frame coding sequences rather
than unprocessed transcripts or genomes. Protein mode uses the original protein
FASTA associated with BUSCO and bypasses CDS preparation.

BUSCO/MetaEuk IDs such as `Species_g123:60-698` map to `Species_g123` in the
original FASTA. The full CDS is translated; coordinates are not used to slice
it because the full table lacks the complete exon/strand model. This reuses the
ortholog assignment but does not reconstruct the exact BUSCO-predicted peptide.
Missing or ambiguous original IDs fail extraction.

## Species sets

```yaml
phylogeny:
  species_sets: [all, phenotyped]
  trait: C4
```

| Set | Species | Output under `results/<analysis>/` |
| --- | --- | --- |
| `all` (default) | All species passing BUSCO and optional species-list selection | `phylogeny/all/` |
| `phenotyped` | Selected species with a nonmissing `phylogeny.trait` | `phylogeny/phenotyped/` |

Traits come from `inputs.species_trait`. Both zero and one are observed; a
single-state or continuous trait can also define the phenotyped subset. Each
inference set needs at least `min_taxa` species (default four). The `all` tree
requires no trait file.

These sets have independent marker plans, alignments, gene trees, and roots.
Coverage is ranked within each set. Changing only `species_sets` requests the
chosen branches without deleting or recomputing unchanged completed branches.
Use a new `analysis` name to retain alternative traits or inference settings.

`phylogeny`, `phylogeny_prepare`, `phylogeny_calibrations`, `timetree`, and
`taxonomy_audit` follow this list. [Molecular-tree contrast pairs](contrast_pairs.md)
do too. The separate `contrast_pairs` target selects and infers a representative
tree under `phylogeny/representatives/`.

## Setup and execution

With Conda deployment, the [phylogeny environment](../workflow/envs/phylogeny.yaml)
supplies cdskit 0.27.0, FAMSA 2.4.1, trimAl 1.5.1, and VeryFastTree 4.0.5.
The [timetree environment](../workflow/envs/timetree.yaml) supplies nwkit 0.27.0
for reference-based rooting and contrast analysis.

Build ASTRAL-IV once with Python 3.12+ and GNU C++ available:

```bash
python workflow/scripts/prepare_phylogeny_tools.py
```

The helper verifies the pinned official source and builds
`resources/phylogeny_tools/bin/astral4_int128` with `LARGE_DATA` (128-bit integers).
It records source/compiler provenance and the executable hash, and reuses a
verified existing build. This build is required above 5,000 species. For offline
compilation, pass `--archives /path/to/archives` containing verified `aster.tar.gz`.

Add the relevant `phylogeny` settings to your dataset configuration, then inspect
the plan or run inference:

```bash
# Resolve the outgroup and select markers, without inferring trees.
./run_pipeline.sh --software-deployment-method conda \
  --configfile config/mydata.yaml --cores 4 --resources mem_gb=16 -- phylogeny_prepare

# Infer gene trees and species trees inside one allocation.
mkdir -p logs
sbatch --partition=YOUR_PARTITION run_pipeline.sh \
  --configfile config/mydata.yaml -- phylogeny
```

The explicit target works with `phylogeny.enabled: false`. Set it to `true` to
include the selected trees in `all`. Stages resume independently; current
manifests determine membership even if files for old markers remain on disk.

## Rooting

`phylogeny.outgroup` accepts `auto` (default) or one exact species label present
in every requested inference set. Automatic rooting selects within each set;
representative selection therefore precedes root selection for `contrast_pairs`.
No extra species is added, and different sets can select different outgroups.

The workflow builds an NCBI guide from the frozen local taxonomy using nwkit.
A binary guide root with a singleton side supplies the outgroup directly. For
an unresolved root, it selects the highest-BUSCO species from each basal lineage
(ties by species ID) and consults nwkit's bundled APG IV order tree. Every basal
lineage must resolve to an order, and the reference must identify a singleton
side that was also a singleton lineage in the input guide.

APG IV is an angiosperm reference. An unresolved root or two multi-species root
sides requires an explicit outgroup. Rooting makes no TimeTree/OpenTree request.
`rooting/outgroup.json` records the candidate manifest, reference checksums, and
root split. The dataset-wide guide lives in `phylogeny/all/rooting/`; phenotyped
inference builds its own guide, and representatives use their compressed guide.

The root is supplied before CASTLES-II length estimation. Its biological
interpretation requires review; the molecular topology itself is not constrained
to NCBI or APG IV.

## Marker and sequence selection

Only unambiguous single-copy `Complete` BUSCO hits are eligible. Duplicated,
fragmented, missing, multiply reported complete hits, and original genes assigned
to multiple markers are omitted. Coverage counts distinct selected species,
so replicate expression runs do not increase a species' weight.

| Setting under `phylogeny` | Default | Use |
| --- | --- | --- |
| `max_markers` | `500` | Cap after ranking eligible markers by coverage descending, then BUSCO ID ascending |
| `min_taxa` | `4` | Minimum species at marker selection and after QC; must be at least four |
| `min_protein_length` | `100` | Known amino acids required per extracted and trimmed sequence |
| `max_unknown_fraction` | `0.05` | Maximum unknown fraction in prepared proteins |
| `trimal_mode` | `gappyout` | `gappyout` or `automated1` |
| `seed` | `12345` | Seed for inference and representative selection |

Coverage is measured before translation and alignment QC, without a minimum
coverage fraction or order-specific condition. Mean BUSCO match length is
recorded but does not affect selection. Loci lost after selection are not
replaced with lower-ranked markers.

CDS preparation calls the pinned cdskit `pad`, `mask`, and `translate` functions,
passing `translation.table` throughout. Module hashes and sequence changes are
recorded. Padding adds N to complete codons and can add bases at the 5′ end to
reduce internal stops. It can change the reading frame and is not proof of a
correct ORF. Original bases are not deleted. U becomes T, X becomes N, and `.`
becomes `-`; other invalid nucleotide symbols fail extraction.

Masking turns stops and unresolved codons into NNN, translated as X, including
terminal stops. Resolvable ambiguity such as GCN → A is retained. Partial-gap
codons are masked; complete-gap codons are translated as gaps and normalized
to X before alignment. Supplied proteins instead have one terminal stop removed,
reject internal stops, and normalize nonstandard residues to X. Known-residue
and unknown-fraction limits then apply.

Single-copy status and sequence QC do not rule out hidden paralogy, gene fusions,
or incorrect translations. Frame changes and masked positions are retained
for inspection in `species/*.json`.

## Alignment and tree inference

FAMSA saves raw alignments in `alignments/raw/`. trimAl selects columns using
`gappyout` by default. `automated1` is available but computes all sequence-pair
identities first; the [method comparison](notes/phylogeny_comparison.md) explains
this choice and its limits.

For column selection, X is temporarily represented as a gap. The saved column
map is applied to the original FAMSA alignment, preserving its residues and X.
Headers, IDs, residues, and column correspondence are checked. Post-trimming
QC removes sequences with fewer than `min_protein_length` known residues and
columns with no known residues among the remaining species.

Markers need `min_taxa` species and at least one variable site with two distinct
standard amino acids. Alignment length and parsimony-informative sites remain
diagnostics; variable loci with zero informative sites can pass. X and gaps
are not observed amino-acid states. Final-to-raw column indices in
`alignments/*.columns.tsv` are 1-based.

Changing only trimAl mode reuses raw FAMSA alignments. Changing the protein-length
threshold also affects extraction and its downstream alignments. Marker ranking
is independent of that length threshold.

VeryFastTree uses double precision and `-lg -gamma`: topology search uses LG+CAT,
then Gamma20 rescales lengths and evaluates likelihoods. This is not a full
LG+Gamma topology search. SH-like local supports are saved, without bootstrap
replicates or support-based filtering.

Gene-tree merging records actual retained coverage and each species' gene-tree
count. ASTRAL requires every selected species in at least one retained gene
tree; zero coverage stops inference and preserves diagnostics. Presence in one
tree is a completeness check, not evidence of adequate phylogenetic information.

ASTRAL-IV estimates topology, local posterior probabilities, and integrated
CASTLES-II substitution lengths using the root, fixed seed, and mean retained
alignment length. No supermatrix is constructed. Missing data, gene-tree error,
hidden paralogy, and model assumptions can affect both topology and lengths.

## Resources

Threads and decimal-GB memory below are per-job scheduling reservations:

| Stage | Thread setting (default) | Memory setting (default GB) |
| --- | --- | --- |
| Preparation | 1 | `preparation_mem_gb: 4` |
| FAMSA | `align_threads: 4` | `alignment_mem_gb: 8` |
| trimAl/QC | 1 | `trimming_mem_gb: 4` |
| VeryFastTree | `tree_threads: 4` | `tree_mem_gb: 8` |
| ASTRAL-IV | `astral_threads: 32` | `astral_mem_gb: 64` |

All settings belong under `phylogeny`. Concurrency follows the total
[resource budget](running.md#resource-budgets). Extraction reads sequences once
per species, then transposes them into per-marker inputs. Measure real-data
benchmarks before scaling up.

## Outputs

Each `results/<analysis>/phylogeny/<set>/` directory contains:

| Output | Contents |
| --- | --- |
| `selection/` (phenotyped only) | Selected sample manifest and trait/source record |
| `rooting/` | Guide tree where applicable, outgroup, and selection evidence |
| `plan/marker_stats.tsv`, `markers.tsv` | Candidate statistics/ranks and selected marker order |
| `plan/species.tsv`, `provenance.json` | Source paths, full-table hashes, lineage, and settings |
| `species/*.faa`, `*.json` | Prepared proteins, sequence QC, cdskit hashes and changes |
| `markers/`, `alignments/raw/` | Per-marker inputs and reusable raw alignments |
| `alignments/*.faa`, `*.json`, `*.columns.tsv` | Trimmed alignments, QC, and final-to-raw columns |
| `gene_trees/`, `gene_trees.nwk`, `gene_trees.json` | Per-marker trees, merged trees, and coverage diagnostics |
| `species_coverage.tsv` | Retained locus counts and representation flags per species |
| `species_tree.nwk`, `species_tree.json` | Rooted tree in substitutions/site and inference provenance |

Logs and resource benchmarks follow the same branch under
`logs/<analysis>/phylogeny/`. Optional [dating](dating.md),
[taxonomic review](taxonomy_audit.md), and [contrast-pair](contrast_pairs.md)
outputs live beside their source tree.

[Automated test coverage](development.md) and [recorded dataset checks](notes/validation.md)
describe what has been validated. The [method comparison](notes/phylogeny_comparison.md)
discusses scientific limitations beyond the processing checks.

## Methods and source documentation

- [GeneGalleon species-tree workflow](https://github.com/kfuku52/genegalleon/blob/main/workflow/core/gg_genome_evolution_core.sh)
- [BUSCO output formats](https://busco.ezlab.org/busco_userguide.html)
- [FAMSA](https://github.com/refresh-bio/FAMSA)
- [cdskit preparation source](https://github.com/kfuku52/cdskit/tree/0.27.0/cdskit)
- [trimAl v1.5.1 source](https://github.com/inab/trimal/tree/d637091abe33595775f40480970d1a18d87a7bcb/source)
- [VeryFastTree](https://github.com/citiususc/veryfasttree)
- [FastTree models, support and Gamma20 scaling](https://morgannprice.github.io/fasttree/)
- [ASTRAL-IV and CASTLES-II, including the >5,000-species build requirement](https://github.com/chaoszhang/ASTER/blob/master/tutorial/astral4.md)
