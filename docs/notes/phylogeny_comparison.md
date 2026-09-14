# GeneGalleon method comparison

[Documentation](../index.md) · [Current phylogeny guide](../phylogeny.md)

The original source review was dated 2026-09-09 and inspected GeneGalleon commit
[`838f237cda4acf384c426a589148e1dbd0029622`](https://github.com/kfuku52/genegalleon/tree/838f237cda4acf384c426a589148e1dbd0029622).
Its container was not run for a performance or accuracy comparison. This note
retains the findings that explain method choices; operational settings and
subsequent implementation changes are documented in the current guides.

## Differences by stage

| Stage | phenoradar_prep | GeneGalleon at the reviewed revision |
| --- | --- | --- |
| BUSCO assignments | Reuses full tables; requires unambiguous single-copy `Complete` hits | Common BUSCO stage creates a summary without the Status column; default selection accepts nonmissing gene IDs without commas |
| Locus selection | Up to 500 eligible loci ranked by species coverage, with deterministic ID ties | No corresponding 500-locus cap; missing/duplicate assignments are removed within each locus |
| CDS preparation | Original CDS with coordinate suffixes removed from IDs; cdskit pad, mask, and translate with an explicit genetic code | Extracts original CDS, uses cdskit pad/mask, then seqkit translation; the inspected pad/mask commands do not explicitly pass a genetic code |
| Alignment | FAMSA protein alignments | MAFFT `--auto`; CDS inputs also have protein-guided backalignment |
| Trimming | trimAl `gappyout`, optionally `automated1`; X treated as missing for selection, original residues restored through column maps | trimAl `automated1`; CDS route uses backtranslation and `-ignorestopcodon` |
| Gene trees | VeryFastTree LG+CAT search with Gamma20 length rescaling | IQ-TREE with default LG+R4; no default per-gene bootstrap option |
| Species topology | ASTRAL-IV with local posterior probabilities | `astral-hybrid --mode 3 --support 2`, subject to the executable distinction below |
| Species-tree lengths | Integrated CASTLES-II, using mean retained locus length | IQ-TREE length optimization on concatenation with a fixed topology; source also contains an unoptimized-tree fallback |
| Dating | LSD2 on the rooted substitution tree, with one estimated rate | CDS route uses IQ2MC and MCMCtree with default IND rate variation; disabled for protein inputs |

Source evidence:
[GeneGalleon defaults](https://github.com/kfuku52/genegalleon/blob/838f237cda4acf384c426a589148e1dbd0029622/workflow/gg_genome_evolution_entrypoint.sh),
[species-tree implementation](https://github.com/kfuku52/genegalleon/blob/838f237cda4acf384c426a589148e1dbd0029622/workflow/core/gg_genome_evolution_core.sh),
[BUSCO table aggregation](https://github.com/kfuku52/genegalleon/blob/838f237cda4acf384c426a589148e1dbd0029622/workflow/support/collect_common_BUSCO_genes.py),
and [sequence extraction](https://github.com/kfuku52/genegalleon/blob/838f237cda4acf384c426a589148e1dbd0029622/workflow/support/batch_extract_busco_fasta.py).
The observation that `collect_common_BUSCO_genes.py` drops Status applies to
this reviewed revision, not all past or future versions.

The ASTRAL comparison depends on the executable installed in GeneGalleon's
environment. Native wASTRAL uses mode 3 for branch-length weighting and mode 4
for unweighted inference. GeneGalleon's `install_astral_hybrid_wrapper` can
instead install an ASTRAL v5 wrapper when `astral-hybrid` is unavailable; that
wrapper discards `--mode` and `--thread`. The command name alone therefore does
not establish weighting or thread behavior. This review did not inspect the
executable inside a running container, so both routes are distinguished.
[wASTRAL modes](https://github.com/chaoszhang/ASTER/blob/master/tutorial/wastral.md),
[GeneGalleon wrapper](https://github.com/kfuku52/genegalleon/blob/838f237cda4acf384c426a589148e1dbd0029622/container/scripts/install_nonconda_fallbacks.sh).

## Sequence and locus assumptions

Both routes translate the full original CDS rather than reconstructing the exact
BUSCO-predicted peptide from its MetaEuk model. The workflow checks ID/lineage
compatibility but lacks checksums proving that the current sequences were the
inputs to the original BUSCO run. Incorrect ORFs, fusions, or extra domains can
remain. BUSCO peptide/domain comparisons would help assess these cases.

cdskit padding can add N at the 5′ end to reduce internal stops, changing the
reading frame. A sequence can consequently pass unknown-residue QC while being
biologically incorrect. The workflow records these changes and tests API results
against the pinned CLI; this verifies implementation rather than ORF identity.
[cdskit pad](https://github.com/kfuku52/cdskit/blob/0.27.0/cdskit/pad.py),
[mask](https://github.com/kfuku52/cdskit/blob/0.27.0/cdskit/mask.py),
[translate](https://github.com/kfuku52/cdskit/blob/0.27.0/cdskit/translate.py).

Coverage ranking favors well-represented loci but does not guarantee shared
information across small clades. The 500-locus cap is a resource choice, not a
demonstrated optimum. Known-residue limits, variable-site checks, and retained
locus counts expose some missingness; presence in one gene tree cannot establish
adequate support for a species placement. Current selection/QC rules are in the
[inference guide](../phylogeny.md#marker-and-sequence-selection).

Single-copy status is assigned within each sample. Ancient duplication and
differential loss can leave different paralogs as apparent single copies across
angiosperms. Hybridization and recombination can also produce histories that
one bifurcating tree does not represent fully.

## Alignment and inference tradeoffs

trimAl `automated1` calculates all sequence-pair identities before choosing
between gappyout and strict. In the pinned source this costs O(N²L) and an N×N
float matrix. At 6,000 species the matrix alone is about 144 MB. Avoiding repeated
all-pairs scans across loci is the main reason for the gappyout default; its
principal gap-statistics scan is O(NL). Equivalent alignment accuracy is not
assumed, and X-as-gap selection adds another difference from the reviewed route.
[trimAl mode selection](https://github.com/inab/trimal/blob/d637091abe33595775f40480970d1a18d87a7bcb/source/autAlignment.cpp),
[column selection](https://github.com/inab/trimal/blob/d637091abe33595775f40480970d1a18d87a7bcb/source/alignment.cpp).

FAMSA supports large protein-family alignment, while VeryFastTree accelerates
FastTree-style approximate search. Their use addresses computational cost but
has not been benchmarked for accuracy on these BUSCO loci. LG+CAT with Gamma20
length rescaling is different from an LG+R4 search. SH-like local support is not
a bootstrap proportion; this workflow does not collapse or weight gene-tree
branches by those support values.
[FAMSA paper](https://www.nature.com/articles/srep33964),
[VeryFastTree paper](https://academic.oup.com/gigascience/article/doi/10.1093/gigascience/giae055/7730000).

ASTRAL-IV/CASTLES-II provides substitution lengths without concatenated length
optimization. Gene-tree estimation error, locus/lineage rate heterogeneity, and
root choice can affect these estimates. Properties established with true gene
trees should be distinguished from performance with finite-sequence estimates.
[ASTRAL-IV documentation](https://github.com/chaoszhang/ASTER/blob/master/tutorial/astral4.md),
[CASTLES paper](https://doi.org/10.1093/bioinformatics/btad221).

## Dating assumptions

TimeTree ages summarize published estimates. Studies can share sequences and
calibrations; a study count is not a count of independent evidence. The adapter
checks sampled-MRCA correspondence and bound consistency, then treats accepted
ranges as hard LSD2 constraints. This differs from fossil priors or MCMCtree
soft tails, and does not verify crown/stem interpretation in each source study.
[TimeTree FAQ](https://timetree.org/faqs),
[nwkit calibration documentation](https://github.com/kfuku52/nwkit/wiki/nwkit-mcmctree).

The LSD2 wrapper estimates one rate, without rate partitions or confidence
simulations. `numsites`, by default the sum of retained trimmed locus lengths,
affects its variance offset; it is not an established effective sample size
for CASTLES lengths. Bounded calibrations can leave multiple equally optimal
absolute scales. Recorded feasible intervals are not statistical confidence
intervals. See [dating](../dating.md) for fitting and output checks, and the
[dated tool evaluation](dating_evaluation.md) for synthetic timings and the
installation rationale.

## Remaining evaluation needs

The workflow would benefit from calibration/rate sensitivity analysis,
comparison of reframed CDS with BUSCO peptides, and branch-specific shared-locus
and quartet diagnostics. Further checks of hidden paralogy, locus independence,
and low-support gene-tree treatment could address errors beyond current QC.

Full-species runs with a small real locus set are needed to measure resources
and information content before broadening inference. The
[recorded dataset checks](validation.md) cover small input/integration examples;
they do not establish full-scale biological accuracy or downstream model validity.
