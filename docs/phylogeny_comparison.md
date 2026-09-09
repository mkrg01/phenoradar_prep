# Phylogeny method comparison with GeneGalleon

Review date: 2026-09-09. This review compares the species-tree workflow in
GeneGalleon commit
[`838f237cda4acf384c426a589148e1dbd0029622`](https://github.com/kfuku52/genegalleon/tree/838f237cda4acf384c426a589148e1dbd0029622)
with the implementation added to phenoradar_prep. GeneGalleon was inspected at
the source level; its container was not run for a performance or accuracy
comparison. Observed code behavior and methodological assessments are
distinguished below.

The implementation provides a reasonable route for reducing the computational
cost of inference for approximately 6,000 species. Gene-tree error, hidden
paralogy, rate heterogeneity and calibration uncertainty remain unresolved
limitations. It has not been demonstrated to be the most accurate method, and
the dated trees should currently be treated as exploratory estimates.

## Comparison by stage

| Stage | phenoradar_prep implementation | Observed GeneGalleon behavior | Assessment and remaining limitations |
| --- | --- | --- | --- |
| BUSCO inputs | Reuses existing full tables and original CDS/proteins; checks lineage name, marker universe and format | Runs BUSCO from CDS/proteins in the common BUSCO stage; reuses results when settings and output records agree | Avoids redundant computation. phenoradar does not have checksums of the inputs used in the original BUSCO run, so it cannot prove that the current CDS files are identical to those inputs |
| Single-copy assignments | Requires one unambiguous `Complete` hit per species and BUSCO; excludes `Duplicated`, `Fragmented`, and original genes assigned to multiple BUSCOs | Defaults to `strictly_single_copy_only=0`; accepts nonmissing gene IDs without commas from a summary that does not retain the BUSCO Status column | Restricting to Complete hits is explicit and conservative. A single Fragmented hit can enter the reviewed GeneGalleon route. Neither rule guarantees orthology |
| Marker selection | Selects **up to 500 loci in descending overall coverage**, among loci with admissible hits in at least four species; ties use lexicographic BUSCO ID order. No coverage floor or order-specific condition; mean BUSCO length is diagnostic | Does not require markers to be present in all species by default; removes species with missing or duplicate assignments from each locus. No 500-locus cap | The rule is reproducible and independent of taxonomic-rank thresholds. The cap is a resource setting, not a demonstrated optimum. Coverage ranking may favor conserved genes and reduce information for shallow divergences |
| Sequence extraction | Maps BUSCO IDs to original CDS IDs and translates the full original CDS, without cropping by MetaEuk coordinates | Removes BUSCO coordinate suffixes and extracts original sequences in batches; uses `cdskit pad` for CDS | Neither route reconstructs the exact BUSCO-predicted peptide. Fusions, incorrect ORFs and extra domains can remain. Corresponding BUSCO peptide outputs or domain checks would improve validation |
| CDS padding and masking | Uses cdskit 0.29.2 `pad` then `mask`; records head/tail N padding, frame changes, masked codon positions and sequence hashes | Uses `cdskit pad` then `cdskit mask`; the inspected commands do not explicitly pass a genetic code to these stages | Avoids terminal-base truncation and rejection of entire CDS sequences for internal stops. Padding can change the frame to reduce stops and does not establish the correct ORF. A mistranslated sequence can still contain few X residues and pass QC |
| Translation and protein QC | Calls cdskit's `translate` implementation and passes `translation.table` to pad, mask and translate. Requires at least 100 known residues and at most 5% unknown residues; preserves masked stops and unresolved codons as X | Translation uses `seqkit translate --allow-unknown-codon --transl-table ...`; backalignment after protein alignment also preserves a DNA analysis route | The processing roles are similar, but the translation programs differ. Tests compare the cdskit API against its CLI. Resolvable ambiguous codons such as GCN are retained. A shared genetic code is used throughout; species- or organelle-specific codes are not supported |
| Alignment | FAMSA 2.4.1, four threads per locus by default | MAFFT `--auto`, one thread per locus; aligns proteins and backaligns CDS inputs to codon sequences | FAMSA is a reasonable choice for large protein alignments. This does not establish that either aligner is uniformly more accurate, and accuracy on these BUSCO loci has not been benchmarked |
| Column trimming | trimAl v1.5.1 `-gappyout` by default, with optional `automated1`; treats X as missing for selection and restores original residues using column maps. No custom 50% gap cutoff | trimAl `-automated1`; the CDS route uses protein-guided backtranslation and `-ignorestopcodon` | gappyout derives a threshold from the gap distribution and avoids automated1's all-pairs identity calculation. It loses the option to switch to strict mode using residue similarity. Treating X as missing is an additional difference, so retained columns need not match |
| QC after trimming | Removes species with fewer than 100 known residues, then all-missing columns. Requires at least four species and observed variation. Records alignment length, variable/informative-site counts and correspondence to raw FAMSA columns | Uses trimAl filtering and subsequent sequence-count checks; not the same thresholds or column-map validation | Separate 100-site and ten-informative-site cutoffs were removed. Variable loci below those counts can remain, without a claim that their trees become more accurate. Nonhomologous columns and partial misalignments can still pass |
| Coverage after QC | Records actual gene-tree coverage and per-species gene-tree counts/distributions without a coverage cutoff. Stops ASTRAL if any selected species is absent from all retained trees | Checks inputs such as a minimum of four gene-tree tips; the species-tree route does not produce the same coverage summary | Measures missingness not visible in the original BUSCO tables. Presence in one tree prevents silent species omission but does not establish adequate information. Joint representation of the four groups around important branches still needs assessment |
| Gene trees | VeryFastTree 4.0.5 with `-double-precision -lg -gamma` and a fixed seed | IQ-TREE with default `LG+R4` and a fixed seed; no default bootstrap option for individual gene trees | VFT searches under LG+CAT and rescales lengths with Gamma20. This differs from IQ-TREE's LG+R4. IQ-TREE offers more extensive search and model options at greater cost for thousands of species and many loci. GeneGalleon does not use ModelFinder by default here |
| Gene-tree uncertainty | Saves SH-like local support but does not collapse low-support branches or use it for weighting | No default individual-gene bootstrap option; species-tree mode 3 requests branch-length weighting. Separate concatenated analyses use 1,000 UFBoot replicates and BNNI | Native wASTRAL can use branch-length weights, unlike the current phenoradar route. SH-like values are not bootstrap proportions. Support definitions and effects on CASTLES should be evaluated before applying common collapse thresholds |
| Species-tree topology | ASTRAL-IV with 128-bit integers; unweighted topology inference and local posterior probabilities | `astral-hybrid --mode 3 --support 2`; native wASTRAL uses branch-length weights, but a bundled ASTRAL v5 compatibility wrapper can change behavior. Separately saves pp1-3, q1-3 and f1-3 annotations | The change can affect weighting as well as speed. Theoretical properties under the multispecies coalescent do not remove errors from estimation, duplication or hybridization. phenoradar currently provides fewer support diagnostics |
| Rooting | Requires a single named outgroup; the current angiosperm configuration uses Amborella and supplies the root before CASTLES | Defaults to taxonomy-based rooting; also supports outgroup, midpoint, MAD and MV, and can transfer the root after length estimation | The selected root needs biological justification and affects CASTLES lengths. The single outgroup's sequence quality and missing data also need assessment |
| Species-tree branch lengths | Integrated CASTLES-II estimates substitutions/site from gene trees, using mean retained alignment length | Fixes the ASTRAL topology and optimizes lengths on concatenated alignments with IQ-TREE `-te ... -n 0`; source also contains a fallback to the unoptimized tree | CASTLES provides substitution units without concatenated ML optimization. VFT length errors and locus/lineage rate heterogeneity can affect the estimates. Unoptimized ASTRAL lengths must not simply be interpreted as substitutions/site |
| Concatenation | Not performed | Estimates concatenated protein/DNA trees and uses concatenation for length optimization after ASTRAL | Omitting concatenation is reasonable under the resource constraints. A concatenated tree can provide a diagnostic under different assumptions, but is not an obligatory reference truth |
| Calibration retrieval | Uses nwkit's HTTP client and TimeTree JSON. Up to 64 representatives and 32 queries, prioritizing larger clades with representatives on every child lineage. No minimum clade size; saves MRCA IDs, studies, missing taxa and raw responses | `nwkit mcmctree --timetree ci --min-clade-prop 0.2`, or manual calibrations | phenoradar retains detailed retrieval and mapping records. Missing representatives in TimeTree can prevent calibration. Prioritizing large clades does not ensure calibration of small clades. nwkit's min-clade-prop filters after retrieval and does not itself reduce query count |
| Dating | CASTLES lengths and calibration bounds feed treePL, with branch-specific rates and an additive smoothing penalty. Runs prime, three random-subsample CV replicates and three final optimizations; no final age CIs | The CDS route uses IQ2MC to construct a likelihood-approximation Hessian from concatenated CDS, then MCMCtree with default IND rate variation. This dating route is disabled for protein inputs | Allowing lineage-rate differences improves applicability to broad angiosperm sampling. Estimates still depend on smoothing, calibration placement and length errors; treePL does not estimate an MCMCtree posterior. Restart age ranges are not CIs. MCMCtree also requires prior and convergence assessment |
| Execution and resources | Snakemake resumes by locus and stage; records memory budgets, benchmarks, inputs and commands. Reuses FAMSA after trimAl setting changes; calls cdskit functions without per-sequence subprocesses | Containers, scheduler support, memory-aware concurrency, artifact provenance and stage-specific ZIP archives | The 8 GB alignment/gene-tree, 4 GB trimAl and 64 GB ASTRAL budgets are scheduling reservations, not measured upper bounds. Full-scale measurements remain necessary; many small files and repeated FASTA reads also have costs |

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

## Rationale and limits

Marker coverage is **the number of species with an admissible single-copy
Complete hit divided by the number of selected unique species**. It is measured
from BUSCO full tables before translation and alignment; replicate runs do not
increase a species' weight. Candidates with hits in at least four species are
ranked by coverage descending, with lexicographic BUSCO ID ties, and the first
500 are selected. All candidates are selected when fewer than 500 are available.
Mean BUSCO match length is diagnostic only. The separate minimum of 100 known
residues applies to individual translated and trimmed sequences. Selection
ranks are saved in `plan/markers.tsv` and `plan/marker_stats.tsv`. Loci lost to
subsequent QC are not replaced, so the final set can contain fewer than 500 loci.

Order-specific species-count and coverage conditions were removed from both
selection and QC. Such thresholds depend on taxonomic-rank definitions and
sampling density. Overall coverage alone still does not guarantee representation
of small clades or shared loci across clades; branch-specific information
diagnostics remain necessary.

The minimum overall coverage of 0.5 was also removed from selection and QC.
Coverage ranking and the locus cap control selection, and coverage after QC
remains recorded. The former post-QC filter could remove loci after species
loss, but the 0.5 boundary had not been validated. Retaining lower-coverage loci
does not remove the need to examine species representation and shared loci.

The restriction of TimeTree queries to clades with at least 20 sampled species
was removed. That count reflects sampling in the input tree, rather than the
number of supporting studies or the reliability of their ages. Representatives
must occur on every child lineage, and workload is bounded by 64 representatives
and 32 queries. Small clades are eligible, but larger clades retain priority.
MRCA correspondence, study counts and age-bound validation remain required.

The separate 100-site alignment cutoff was redundant with the default minimum
of 100 known residues per retained sequence. The ten-informative-site and
ten-gene-trees-per-species cutoffs were also removed because they had not been
validated for these data. Informative-site counts remain in alignment QC JSON;
per-species gene-tree counts and their distribution remain in
`species_coverage.tsv` and `gene_trees.json`. A variable site contains at least
two observed standard amino-acid states. Variable sites need not be
parsimony-informative, so variable loci with zero informative sites can remain.
Invariant loci are excluded, and any species absent from all retained gene
trees stops inference. These changes reduce arbitrary exclusions without
establishing the accuracy of low-information gene trees. The per-sequence
100-known-residue and 5%-unknown limits and TimeTree's five-study condition
remain in place.

CDS preparation follows GeneGalleon's approach and calls the actual cdskit
`pad`, `mask` and `translate` functions. Sources are pinned, API-boundary hashes
are checked, and tests compare results with the CLI. Original input files remain
unchanged. JSON records include hashes after padding and masking, head/tail
padding counts and masked positions. Masking allows translation to continue
across internal stops, retaining unresolved codons as X. Known-residue counts
and unknown fractions exclude sequences dominated by missing data. Terminal
stops also become X and contribute to column selection. Protein inputs bypass
CDS preparation and are rejected for internal stops.
[cdskit pad](https://github.com/kfuku52/cdskit/blob/218f6ed61abcac7f11dd81b17c087cb16119c39e/cdskit/pad.py),
[mask](https://github.com/kfuku52/cdskit/blob/218f6ed61abcac7f11dd81b17c087cb16119c39e/cdskit/mask.py),
[translate](https://github.com/kfuku52/cdskit/blob/218f6ed61abcac7f11dd81b17c087cb16119c39e/cdskit/translate.py).

Padding does more than make sequence lengths divisible by three. When internal
stops occur, it can add N at the 5' end and choose a frame with fewer stops.
For a correctly predicted CDS containing a real internal stop or frameshift,
this can retain an unrelated translation. Comparing reframed sequences with
BUSCO peptides or domains is a priority; successful padding is not evidence
that the resulting sequence is biologically correct.

trimAl `automated1` first computes identities for all sequence pairs and then
chooses gappyout or strict. The pinned v1.5.1 code allocates an N-by-N float
matrix and performs O(N²L) work in this step. For 6,000 species, the matrix alone
is approximately 144 MB (6,000² times four bytes); this allocation alone is not
the main anticipated memory bottleneck. Avoiding repeated all-pairs scans
across many loci is the principal reason for the gappyout default. Its main
gap-statistics scan is O(NL). It omits the residue-similarity-based switch to
strict, so equivalent accuracy is not assumed.
[automated1 selection](https://github.com/inab/trimal/blob/d637091abe33595775f40480970d1a18d87a7bcb/source/autAlignment.cpp),
[gappyout and column mapping](https://github.com/inab/trimal/blob/d637091abe33595775f40480970d1a18d87a7bcb/source/alignment.cpp).

The temporary trimAl input replaces X with gaps. Its `-colnumbering` output
selects columns from the original alignment, preserving X and known residues.
`-keepseqs` and explicit ID-set validation leave species removal to the QC
stage. The workflow does not use `nogaps`: with thousands of species, missing
data or indels in a few species could otherwise remove many columns. These
procedures do not directly identify misalignments, fused genes or hidden
paralogs.

FAMSA was developed for aligning large protein families, supporting its use
here. Its published family-alignment benchmarks do not directly establish
accuracy for these BUSCO inputs.
[FAMSA implementation](https://github.com/refresh-bio/FAMSA),
[FAMSA paper](https://www.nature.com/articles/srep33964).

VeryFastTree parallelizes and accelerates FastTree-style approximate search.
Specifying LG does not make LG+CAT and LG+R4 searches equivalent. The `-gamma`
option must not be described as performing a full LG+Gamma ML topology search.
[VeryFastTree paper](https://academic.oup.com/gigascience/article/doi/10.1093/gigascience/giae055/7730000).

ASTRAL-IV handles many gene trees with missing species and integrates CASTLES-II
length estimation in substitution units. This provides the units needed for
dating without concatenated length optimization. Theoretical properties given
true gene trees should be distinguished from accuracy with trees estimated
from finite sequences.
[ASTRAL-IV documentation](https://github.com/chaoszhang/ASTER/blob/master/tutorial/astral4.md),
[CASTLES paper](https://doi.org/10.1093/bioinformatics/btad221).

BUSCO single-copy status is assigned within each sample. Ancient whole-genome
duplication and differential loss can leave different paralogs as single copies
in different angiosperms, which Complete-only filtering does not detect.
Hybridization, gene flow and recombination can also produce histories that one
bifurcating species tree cannot fully describe. The current implementation does
not resolve these issues.

TimeTree ages are secondary calibrations derived from published estimates.
Studies can share data and fossil calibrations, so a reported count of 35
studies does not imply 35 independent pieces of evidence. Converting ranges
of variation among studies into hard treePL bounds is not equivalent to using
soft priors in MCMCtree. Matching sampled MRCAs does not validate every source
study's crown/stem interpretation.
[TimeTree FAQ](https://timetree.org/faqs),
[nwkit calibration functionality](https://github.com/kfuku52/nwkit/wiki/nwkit-mcmctree).

treePL uses penalized likelihood with branch-specific rates, replacing a common
rate across the tree. Smoothing is selected by mean scores across repeated
random-subsample CV runs, followed by repeated final optimizations at that
smoothing value. Native LF/PL/CV iteration limits bound work. The `thorough`
option is not the default because small tests repeatedly returned to the same
solution while continuing to the iteration limit. Objective/age reproducibility
and native convergence flags are recorded separately; reproducibility does not
prove global convergence.
[treePL paper](https://doi.org/10.1093/bioinformatics/bts492),
[treePL run options](https://github.com/blackrim/treePL/wiki/Run-Options).

A 2026 study reported CASTLES-Pro plus treePL length estimation and dating for
10,000 species and 1,000 genes, averaging 1.6 hours and at most 26 GB for those
stages. These figures exclude gene-tree and species-tree topology inference.
The study's conditions differ from this CASTLES-II protein pipeline. Following
that study, the default treePL `numsites` is the sum of retained, trimmed gene
lengths; this does not establish the effective sample size of CASTLES lengths
with missing data and variation in locus rates.
[Tabatabaee et al. (2026)](https://doi.org/10.1093/sysbio/syag038).

The pinned treePL source contains defects including an uninitialized random-CV
index, an uninitialized prime flag and inconsistent CV optimizer selection.
The build applies documented corrections. Methodological justification and
implementation correctness require separate checks.
[Patch details and rationale](../workflow/patches/README.md).
The workflow also records branches that treePL raises to `1/numsites`.

## Example Methods text

This example describes the defaults. Supply the actual species and retained
locus counts, genetic code and software versions from the run records, and
follow it with descriptions of tree inference and dating matching the stages
and settings above.

Single-copy Complete BUSCO hits were identified from the existing full tables
using a common lineage dataset. Ambiguous assignments, including original genes
assigned to multiple BUSCO markers, were excluded. Loci present in at least four
selected species were ranked by species coverage, with
ties resolved lexicographically by BUSCO identifier. Up to 500 loci were selected;
loci excluded subsequently were not replaced. Mean BUSCO match length was
recorded without a selection cutoff. Original CDS sequences were processed using
cdskit pad, mask and translate with the same specified genetic code throughout.
Proteins with fewer than 100 known amino acids or more than 5% unknown residues
were excluded. Loci were aligned using FAMSA and trimmed using trimAl gappyout.
For column selection, X was treated as a gap; selected columns were then extracted
from the original alignment to preserve residue identities. After trimming,
sequences with fewer than 100 known residues and columns lacking known residues
were removed. Loci with fewer than four remaining species or no observed variable
sites were excluded. Alignment length and parsimony-informative-site counts were
recorded without additional numerical cutoffs. Species coverage of retained gene
trees was recorded without a minimum coverage fraction. Per-species gene-tree counts were
recorded, and species-tree inference required every selected species to be
represented in at least one retained gene tree.

For automatic TimeTree retrieval, also describe its selection rules:

Up to 64 representative species were selected by repeatedly partitioning the
largest clade and choosing the species with the highest gene-tree count in each
part. Internal nodes with representatives on every child lineage were eligible
for TimeTree queries, without a minimum clade size. Up to 32 eligible nodes were
queried in descending order of clade size, with ties resolved lexicographically
by representative species labels. Returned taxa had to represent every child
lineage and identify the intended MRCA. Bounds supported by at least five
distinct bibliographic records were required to be finite and positive, include
the reported age, and be consistent with the other retained calibrations.

## Improvement priorities

1. **Assess calibration and optimization sensitivity.** treePL supports rate
   variation. Expand the CV grid when its optimum is at an endpoint, and examine
   objective and age differences among restarts. Assess removal or replacement
   of major calibrations, `numsites`, and short-branch adjustments. Propagating
   calibration, gene-tree and branch-length uncertainty remains unimplemented.
2. **Improve sequence and locus diagnostics.** Report unusual lengths, unknown
   residues, anomalous alignment regions and terminal branches. Automatic
   long-branch removal can discard genuinely rapid evolution, so investigate
   causes first. cdskit and trimAl are integrated, but reframed sequences still
   need comparison with BUSCO peptides or the corresponding domains.
3. **Report information supporting species-tree branches.** Include competing
   quartet support, effective gene counts per branch and shared loci across
   clades alongside localPP. This is a higher priority than replacing the
   500-locus cap with another arbitrary count.
4. **Address gene-tree estimation error.** Evaluate low-support branch handling
   and integration using support or branch lengths, including the weighting in
   GeneGalleon's native wASTRAL route. Do not interpret SH-like values as
   bootstrap proportions. More extensive searches for selected loci could
   control cost, but selection must not favor a desired topology.
5. **Assess independence and hidden paralogy.** Genomic locations and gene-family
   information from representative genomes could help. Removing multiple
   assignments of the same original gene does not assess linkage or ancient
   duplication.
6. **Measure resources at the full species count.** Use a small set of real loci
   across all species to measure RAM, runtime, alignment lengths, missingness
   and anomalous branches. This validates the selected pipeline's feasibility
   without requiring a comparison of numerous software packages.

## Implemented changes and validation status

The BUSCO phylogeny branch integrates cdskit padding, masking and translation
with trimAl. It replaces terminal-base truncation and whole-CDS rejection for
internal stops, shares the genetic code across all preparation steps, and
records frame changes, masked positions and column correspondence. FAMSA and
trimAl are separate restartable stages. Marker selection uses overall coverage
ranking, deterministic ID ties and a maximum of 500 loci. It does not use
order-level conditions or length-based ranking.

Mean BUSCO length, alignment length, informative-site counts and per-species
gene-tree counts are retained as diagnostics after removal of the former
100-aa, 100-site, ten-site and ten-gene selection cutoffs. The separate
per-sequence known-residue minimum remains explicit. cdskit's padding,
masking and translation behavior is preserved, including possible 5' padding;
no additional reading-frame restriction is imposed.

TimeTree retrieval through nwkit includes caching, taxon matching and age-bound
validation and feeds treePL dating. treePL is the sole dating engine; the
earlier engine and rate-partition configuration were removed. Minimum overall
coverage and TimeTree clade-size cutoffs were removed while their diagnostic
values remain recorded. ASTRAL's 128-bit build provenance and executable hash
are checked, so renaming an ordinary binary `astral4_int128` cannot bypass the
requirement for more than 5,000 species.

Tests cover sequence/tree processing contracts, incorrect calibration mappings,
missing taxa, repeated MRCAs, ancestor-age conflicts, HTTP failures, offline
reuse, cdskit CLI agreement, and a small workflow using actual FAMSA, trimAl,
VeryFastTree, ASTRAL and treePL. They also check FAMSA reuse after trimAl setting
changes, retention of low-coverage loci, and small-clade calibration eligibility
under a bounded query budget.

A real-input check used 20 shared Complete BUSCOs from tlight for Amborella
trichopoda, Oryza sativa, Abelia chinensis and Abeliophyllum distichum. All 80
CDS passed translation QC, without a frame change from head padding; terminal
stops were masked. The previous QC retained 18 loci and excluded two with eight
informative sites each. Repeating trimAl/QC on the same raw alignments with the
simplified conditions retained all 20. The recovered loci have alignment lengths
of 189/150 and variable-site counts of 134/113, respectively. The original 18
retained alignments and column maps were unchanged. This checked input
compatibility through trimming; it did not infer species or dated trees from
these four species. Original and updated summaries are stored locally in
`results/phylogeny_cdskit_trimal_smoke/summary.json` and
`results/phylogeny_qc_simplification_smoke/summary.json`. These 20 loci had been
preselected for match lengths of at least 100 aa, so this check did not assess
all candidates with shorter matches.

A separate live TimeTree check made three queries for a four-species test tree
and retrieved calibration candidates and study information. Its branch lengths
were artificial, not inferred from BUSCO.

Inference has not been run for all 5,586 species. Full-scale memory requirements,
biological accuracy, calibration sensitivity, node-age uncertainty and effects
on downstream conclusions such as PGLS remain untested. Successful software
tests do not establish biological validity.
