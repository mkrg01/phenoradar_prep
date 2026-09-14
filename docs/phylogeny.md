# BUSCO species trees

[Documentation](index.md) · [Dating](dating.md)

The `phylogeny` target infers a rooted species tree from existing BUSCO full
tables and original CDS or proteins, independently of ODB/KEGG analyses.

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

Use the original sequences from each BUSCO run. CDS must be oriented, in-frame
coding sequences rather than unprocessed transcripts or genomes. Protein mode
uses the corresponding protein FASTA and bypasses CDS preparation.

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

| Set | Species | Output under `results/<run_name>/` |
| --- | --- | --- |
| `all` (default) | All species passing BUSCO and optional species-list selection | `phylogeny/all/` |
| `phenotyped` | Selected species with a nonmissing `phylogeny.trait` | `phylogeny/phenotyped/` |

Traits come from `inputs.species_trait`. Both zero and one are observed; a
single-state or continuous trait can also define the phenotyped subset. Each
inference set needs at least `min_taxa` species (default four). The `all` tree
requires no trait file.

Each set has independent markers, alignments, trees, and roots. All phylogeny,
dating, audit, and molecular-pair targets follow `species_sets`. The separate
`contrast_pairs` target infers a representative tree. Use a new `run_name` to
retain alternative traits or settings.

## Setup and execution

Analysis tools are bundled in the container. Native execution uses
[phylogeny.yaml](../workflow/envs/phylogeny.yaml) and
[timetree.yaml](../workflow/envs/timetree.yaml); ASTRAL-IV is prepared automatically.

```bash
# Select markers and resolve the outgroup.
./run_pipeline.sh --configfile config/mydata.yaml \
  --cores 4 --resources mem_gb=16 -- phylogeny_prepare

# Infer gene trees and species trees.
./run_pipeline.sh --configfile config/mydata.yaml \
  --cores 32 --resources mem_gb=128 -- phylogeny
```

Set `phylogeny.enabled: true` to include trees in `all`. For offline native ASTRAL
setup, run `python workflow/scripts/prepare_phylogeny_tools.py --archives /path/to/archives`
with the verified `aster.tar.gz`, Python 3.12+, and GNU C++. The helper builds the
pinned 128-bit ASTRAL-IV executable required for datasets above 5,000 species.

## Rooting

`phylogeny.outgroup` accepts `auto` or one exact species ID present in every
requested set. Automatic selection uses the local NCBI guide, with nwkit's APG IV
order tree as a fallback for angiosperms. It chooses within each set and adds no
species. If no unambiguous outgroup is found, supply one explicitly.

Review `rooting/outgroup.json` for the choice and supporting guide. The outgroup
is supplied before CASTLES-II length estimation; molecular topology is not
constrained to the guide.

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

Coverage is ranked before sequence/alignment QC; loci lost later are not replaced.
There is no minimum coverage fraction or order-specific condition.

CDS preparation uses cdskit `pad`, `mask`, and `translate` with `translation.table`.
Padding can change the reading frame; stops and unresolved codons become X.
Supplied proteins have one terminal stop removed, reject internal stops, and
normalize nonstandard residues to X. Inspect changes and masking in `species/*.json`:
single-copy status and QC do not establish a correct ORF or exclude hidden paralogy.

## Alignment and tree inference

FAMSA saves raw alignments, then trimAl selects columns (`gappyout` by default;
`automated1` is available). X is treated as a gap for column selection, then the
selected columns are recovered from the original alignment, preserving X.
Final-to-raw column indices in `alignments/*.columns.tsv` are 1-based.

After trimming, sequences need `min_protein_length` known residues. Columns with
no known residues are removed. Each locus needs `min_taxa` species and at least
one variable amino-acid site; parsimony-informative sites are diagnostic only.
Changing trimAl mode reuses raw FAMSA alignments.

VeryFastTree uses double precision and `-lg -gamma`: LG+CAT topology search with
Gamma20 length rescaling. SH-like local supports are saved without bootstrap
replicates or support filtering.

ASTRAL-IV combines gene trees and estimates local posterior probabilities and
CASTLES-II branch lengths in substitutions/site. Every selected species must
occur in at least one retained gene tree; inspect `species_coverage.tsv` for
uneven coverage. No supermatrix is constructed. Missing data, gene-tree error,
paralogy, and model assumptions can affect topology and lengths.

## Resources

Threads and decimal-GB memory below are per-job scheduling reservations:

| Stage | Thread setting (default) | Memory setting (default GB) |
| --- | --- | --- |
| Preparation | 1 | `preparation_mem_gb: 4` |
| FAMSA | `align_threads: 4` | `alignment_mem_gb: 8` |
| trimAl/QC | 1 | `trimming_mem_gb: 4` |
| VeryFastTree | `tree_threads: 4` | `tree_mem_gb: 8` |
| ASTRAL-IV | `astral_threads: 32` | `astral_mem_gb: 64` |

All settings belong under `phylogeny`; see [resource budgets](running.md#resource-budgets)
for concurrency and allocation sizing.

## Outputs

Each `results/<run_name>/phylogeny/<set>/` directory contains:

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
`logs/<run_name>/phylogeny/`. Optional [dating](dating.md),
[taxonomic review](taxonomy_audit.md), and [contrast-pair](contrast_pairs.md)
outputs live beside their source tree.

## Methods and source documentation

- [GeneGalleon species-tree workflow](https://github.com/kfuku52/genegalleon/blob/main/workflow/core/gg_genome_evolution_core.sh)
- [BUSCO output formats](https://busco.ezlab.org/busco_userguide.html)
- [FAMSA](https://github.com/refresh-bio/FAMSA)
- [cdskit preparation source](https://github.com/kfuku52/cdskit/tree/0.27.0/cdskit)
- [trimAl v1.5.1 source](https://github.com/inab/trimal/tree/d637091abe33595775f40480970d1a18d87a7bcb/source)
- [VeryFastTree](https://github.com/citiususc/veryfasttree)
- [FastTree models, support and Gamma20 scaling](https://morgannprice.github.io/fasttree/)
- [ASTRAL-IV and CASTLES-II, including the >5,000-species build requirement](https://github.com/chaoszhang/ASTER/blob/master/tutorial/astral4.md)
