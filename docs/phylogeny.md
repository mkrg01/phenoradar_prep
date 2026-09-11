# BUSCO protein phylogeny and optional dating

[Back to README](../README.md)

The optional `phylogeny` target reuses BUSCO results to infer a rooted species
tree. It does not run ODB-mapper, concatenate alignments, optimize branch lengths
on a concatenated alignment, or run IQ-TREE/MCMCtree.

```text
Selected species + existing BUSCO full tables + original in-frame CDS/proteins
  -> select BUSCO markers -> extract CDS once per species
  -> cdskit pad -> mask -> translate (or supplied proteins)
  -> transpose to per-marker FASTA -> FAMSA -> trimAl gappyout -> alignment QC
  -> VeryFastTree -> ASTRAL-IV int128 + integrated CASTLES-II
  -> species_tree.nwk (substitutions/site)
  -> optional manual or nwkit/TimeTree secondary calibrations + LSD2
  -> species_tree.dated.nwk (million years)
```

## Species sets and reusable outputs

`phylogeny.species_sets` selects one or both independent inference runs:

```yaml
phylogeny:
  species_sets: [phenotyped]  # default: [all]; both: [all, phenotyped]
  trait: C4
```

| Set | Inference species | Directory under `results/<analysis>/` |
| --- | --- | --- |
| `all` | All species passing the ordinary BUSCO/species-list selection | `phylogeny/` (unchanged) |
| `phenotyped` | Those selected species with a nonmissing value in `phylogeny.trait` | `phylogeny_phenotyped/` |

Traits come only from `inputs.species_trait`; metadata trait columns are ignored.
The existing name normalization and missing-value rules apply: blanks, `NA`,
`NaN`, `nan`, and absent species rows are unknown. Both `C4=0` and `C4=1` are
observed. All observed species are retained, with no representative skim and
no requirement for two states. The selected set must contain at least
`phylogeny.min_taxa` distinct species (default four). `phylogeny.trait` is
independent of `contrast.trait`. The `all` run needs no trait file.

The `phylogeny`, `phylogeny_prepare`, `phylogeny_calibrations`, and `timetree`
targets all follow `species_sets`. `phylogeny.enabled: true` also includes the
requested sets in `all`; it is unnecessary for the explicit targets.
To try observed species first, run `phylogeny` with `[phenotyped]`. Later change
the list to `[all]` or `[all, phenotyped]` and run the same target. Switching this
list alone does not delete or recompute completed inference or dating outputs.
Each directory retains its own marker plan, alignments, gene trees, rooting,
calibrations, dating, and benchmarks; logs use matching branch directories.
Changing source inputs or inference settings still triggers the affected jobs.
The shared `run.json` records the latest configuration; branch provenance and
`phylogeny_phenotyped/selection/selection.json` retain the relevant input records.
Use different `analysis` names to retain multiple trait definitions or inference
settings simultaneously.

The observed-species run does not depend on a completed all-species tree or
alignment. Marker coverage is ranked separately within each species set, and
alignments and trees are inferred independently using the same rules and settings.
`contrast_pairs` remains a separate representative-selection analysis.

Manual calibration taxa and an optional TimeTree representative list must be
valid for every requested tree. Automatic TimeTree representatives and
calibrations are selected separately for each tree; raw API responses share
the existing cache.

## Input requirements

Set `phylogeny.busco_full_dir` to a directory of per-species BUSCO full tables.
The default filename is `<species>.busco.full.tsv`; `.tsv.gz` is supported when
configured in `busco_full_suffix`. A summary containing only S/D/F/M counts is
not sufficient. Full tables must have `Busco id`, `Status`, `Sequence`, `Score`,
and `Length` columns, a lineage header, and the same complete marker-ID set.
Use one lineage dataset/version consistently. Genome-mode full tables are not
accepted because reconstructing spliced proteins requires more information.

By default, proteins are translated from the original, already oriented,
in-frame CDS in the selected metadata's `samples.tsv`. An optional
`sequence_dir` and `sequence_suffix` select a different original sequence
directory. `sequence_mode: protein` requires `sequence_dir` and a matching
protein filename suffix. Do not supply unprocessed transcripts or genomes as CDS.

For the existing tlight data, the required pair is:

- `tree_all/gfe_data/busco_full_longest_cds/<species>.busco.full.tsv`
- `tree_all/gfe_data/longest_cds/<species>_longestCDS.fa.gz`

BUSCO 5.8.2/MetaEuk transcriptome hits such as `Species_g123:60-698` are
mapped back to `Species_g123`. The full original CDS is translated, just as in
the reference workflow. **The coordinates are not used to slice the CDS**:
full tables do not retain the complete MetaEuk exon/strand model. This reuses
BUSCO's ortholog assignment; it does not reconstruct the exact BUSCO-predicted
peptide. Original CDS IDs must match. Missing or ambiguous IDs fail the job.

Only unambiguous `Complete` (single-copy) hits are retained. Duplicated,
fragmented, missing, multiply assigned original genes, and multiply reported
complete hits are omitted. This reduces obvious ambiguity but does not prove
orthology in the presence of ancient duplication and differential loss.
CDS preparation uses Bioconda cdskit 0.27.0: `pad` -> `mask` -> `translate`.
The workflow calls the same Python functions as those commands to avoid starting
three subprocesses for every CDS. The three module checksums are verified before
use and recorded in each species' JSON. `translation.table` is passed explicitly
to all three operations. This affects the BUSCO phylogeny branch only; the
ODB/KEGG translation rules retain their existing behavior.

`pad` adds `N` to complete codons, and may also add bases at the **5' end** to
reduce internal stops, even when the original length is divisible by three.
This adopts GeneGalleon's heuristic; it does not establish the biological ORF
or correct internal frameshift errors. Head/tail padding, reading-frame changes,
the native padding report, and padded/masked sequence hashes are retained.
No original nucleotides are deleted by padding. Input U becomes T, X becomes N,
and `.` becomes `-`; other invalid nucleotide symbols fail extraction.

`mask` replaces stop codons and ambiguous codons without a resolved amino acid
with `NNN`, which translates to `X`. Resolvable ambiguity such as `GCN` -> A is
retained. Partial-gap codons are masked; complete-gap codons are translated as
gaps, then normalized to X before alignment. Translation continues across masked
codons; terminal stops also become X rather than truncating the sequence. The
audit records changed codon positions and residual internal stops before masking.
Proteins require at least 100 **known** amino acids and at most 5% unknown
residues by default. These checks limit missing data but do not validate a
reframed sequence. Supplied proteins bypass CDS preparation: one terminal `*`
is removed, internal `*` causes rejection, and the same known-residue/unknown
limits apply. Nonstandard amino-acid symbols are normalized to X.

## Setup and execution

Run from the repository root. The launcher with
`--software-deployment-method conda` creates the required stage environments:

- [`phylogeny.yaml`](../workflow/envs/phylogeny.yaml): cdskit 0.27.0, FAMSA
  2.4.1, trimAl 1.5.1, VeryFastTree 4.0.5, and Python dependencies from Conda.
- [`timetree.yaml`](../workflow/envs/timetree.yaml): nwkit 0.27.0 from Conda
  for automatic outgroup selection, contrast skims/figures, and optional calibration retrieval.
- [`dating.yaml`](../workflow/envs/dating.yaml): Conda Python and compiler
  dependencies; its adjacent [post-deploy script](../workflow/envs/dating.post-deploy.sh)
  automatically builds unmodified upstream LSD2 2.4.4 into the environment.

LSD2 uses official commit `c61110f3a4fa05325b45c97b2134792ff9d55d4c`,
a verified source archive, and the upstream makefile with Conda's GNU C++ 13
compiler. It needs no manual installation or project source patch. The build
record is saved as `share/lsd2/build.json` inside the environment. For an offline
source build, set `LSD2_SOURCE_ARCHIVE=/path/to/lsd2.tar.gz` to the archive for
that commit; the SHA-256 check still applies. Conda dependencies must also be
available locally for an entirely offline installation.

IQ-TREE contains LSD2, but its documented dating interface also takes an
alignment. A tree-only test with Bioconda IQ-TREE 3.1.3 did not perform dating.
The standalone executable accepts the existing CASTLES substitution tree
directly. See the [evaluation](dating_evaluation.md) for the tested invocation
and its limitations.

Build ASTRAL-IV once with Python 3.12+ and GNU C++:

```bash
python workflow/scripts/prepare_phylogeny_tools.py
```

The helper downloads the fixed official ASTRAL source, checks its SHA-256,
and builds `resources/phylogeny_tools/bin/astral4_int128`. Build provenance
and executable checksums are recorded next to `bin/`. It follows the official
`LARGE_DATA` (128-bit integers) instruction without machine-specific flags.
This matters above 5,000 species. An existing verified build is reused.
An offline build accepts `--archives /path/to/archives`, containing the verified
`aster.tar.gz`. This helper no longer builds trimAl or treePL.

Create an override, e.g. `config/phylogeny.local.yaml`:

```yaml
phylogeny:
  busco_full_dir: /path/to/tlight/tree_all/gfe_data/busco_full_longest_cds
  sequence_dir: /path/to/tlight/tree_all/gfe_data/longest_cds
  outgroup: auto  # automatic; or an exact selected species label
```

The default `auto` selects an outgroup **within the species set
actually used by each inference run**. The full and phenotyped branches use
their respective input manifests; the contrast branch first selects
representatives with nwkit skim, then chooses among those representatives.
No species is added for rooting, and
the branches can select different outgroups.
Use `auto` for automatic selection; `null` is rejected with a configuration error.
`nwkit constrain` builds an NCBI guide with the existing frozen SQLite snapshot;
contrast rooting uses the compressed guide containing exactly its representatives.
If its root is binary with a singleton side, that species is used. For an
unresolved root, the workflow takes the highest-BUSCO species from each basal
lineage (ties by species ID) and consults nwkit's bundled **APG IV order tree**.
Exact taxids from the snapshot map these representatives to APG IV orders.
The reference must resolve every basal lineage and identify a singleton side
that was also a singleton lineage in that run's NCBI guide. This avoids
mistaking one representative of a multi-species outgroup for an outgroup to
all other species. APG IV is an angiosperm reference, not a universal guide.
Unresolved roots or two multi-species root sides require an explicit outgroup.
No OpenTree/TimeTree request is made for rooting. The root position is a
reference-based choice, not an ancestral-root estimate from BUSCO sequences.

`results/<analysis>/rooting/` contains the dataset-wide `ncbi_tree.nwk` and
`taxids.tsv`. The phenotyped branch builds its own guide from its selected
species under `phylogeny_phenotyped/rooting/`. Each inference branch writes its
own `rooting/outgroup.txt` and `rooting/outgroup.json` under `phylogeny/`,
`phylogeny_phenotyped/`, or `contrast/phylogeny/` (candidate
manifest/count, source, reference checksum and root split). An explicit outgroup
bypasses automatic selection and must already belong to that run's inference
manifest; a name absent from the skim representatives fails rather than adding
a species. All runs resolve their own outgroup **before CASTLES-II length estimation**.
The molecular topology is not constrained to the NCBI/APG IV topology.

```bash
# Check the graph and inputs.
./run_pipeline.sh --configfile config/mydata.yaml config/phylogeny.local.yaml \
  --cores 32 --resources mem_gb=128 --dry-run -- phylogeny

# Select markers and audit all BUSCO inputs without estimating trees.
./run_pipeline.sh --software-deployment-method conda \
  --configfile config/mydata.yaml config/phylogeny.local.yaml \
  --cores 4 --resources mem_gb=16 -- phylogeny_prepare

# Infer gene trees and the species tree within one Slurm allocation.
sbatch --partition=YOUR_PARTITION run_pipeline.sh \
  --configfile config/mydata.yaml config/phylogeny.local.yaml -- phylogeny
```

`phylogeny.enabled: true` adds the branch to the full workflow's `all` target.
The explicit target works while it is false. Input files are read in place.
Stages resume independently, and changed selected species are propagated through
the marker plan. Old files for deselected markers/species are not included in
the current inference; the current manifests define membership.

A step-by-step scientific review against a fixed GeneGalleon revision, including
limitations and priorities, is in [the method comparison](phylogeny_comparison.md).

The optional [`taxonomy_audit` target](taxonomy_audit.md) uses MonoPhy on the existing
species tree and NCBI taxonomy to report non-monophyletic groups, native
intruder/outlier flags, associated run IDs and figures at each rank. It follows `species_sets` and changes no
inference outputs or contrast pairs. Species removal is not part of this stage.

## Selection and inference defaults

The default selects the **500 eligible BUSCO markers with the highest overall
coverage**, or all eligible markers if fewer than 500 are available. Coverage is
the number of selected species with an admissible single-copy `Complete` hit,
divided by the total number of selected unique species. It is measured from the
existing BUSCO full tables, before translation and alignment QC. Multiple runs
for one species do not increase its weight. Duplicated/fragmented/ambiguous hits
do not count toward coverage.

Eligibility requires at least four species (`min_taxa: 4`), without a minimum
coverage fraction. Eligible markers are sorted by coverage descending;
equal coverage is resolved by BUSCO marker ID in ascending lexicographic order.
Mean BUSCO match length is recorded as a diagnostic and does not affect selection.
The separate `min_protein_length: 100` filter applies to individual extracted and
trimmed sequences, not BUSCO match lengths. `max_markers: 500`
sets the cap. These are adjustable starting settings, not a guarantee that 500
loci resolve every branch.

Selection does not use order membership or order-level coverage. The obsolete
settings `clade_rank`, `clade_min_species`, and `min_clade_occupancy` have been
removed; delete them from local overrides if present. `plan/markers.tsv` lists
the selected markers in selection order, with a 1-based `selection_rank`.
`plan/marker_stats.tsv` includes ranks for all eligible markers (blank for
ineligible markers), and `plan/provenance.json` records the ranking rule and
coverage definition. Loci excluded by later QC are not replaced with lower-ranked
markers, so the final gene-tree set can contain fewer than 500 loci.

FAMSA 2.4.1 uses its default guide-tree strategy. Raw alignments are saved under
`alignments/raw/`. trimAl selects columns using `-gappyout` by default: an
adaptive gap-distribution threshold, not removal of every column containing a
gap. `phylogeny.trimal_mode: automated1` is also available. GeneGalleon uses
`automated1`; that mode first computes all sequence-pair identities, taking
O(N²L) work and an N×N float matrix. The default `gappyout` avoids this step,
using gap statistics whose principal scan is O(NL). Neither mode proves that
retained columns are homologous. The former `max_gap_fraction` setting is removed
and rejected if supplied, so an obsolete threshold cannot be silently ignored.

For selection, a temporary copy represents X as `-`, so masking contributes to
gap statistics. The retained column map is then applied to the original FAMSA
alignment, preserving its X and known residues. `-keepheader -keepseqs` plus
explicit ID, residue and column-map validation prevent silent taxon changes.
This missing-data handling differs from directly running GeneGalleon's trimAl
command, even when `automated1` is selected. Post-trimming QC removes taxa with
fewer than 100 known residues and columns with no known residues among the
remaining taxa. No additional 50% gap filter is applied. Markers require four
species and at least one observed variable site (two or more distinct standard
amino acids); invariant/all-missing loci are excluded. Alignment length,
variable sites and parsimony-informative sites are recorded in the alignment QC
JSON. There is no separate minimum alignment length or informative-site cutoff.
Variable loci with zero parsimony-informative sites are allowed; their retention
does not imply that their topology is well resolved. X and gaps do not count as
observed amino-acid states. All dropped markers and species are recorded.
`alignments/*.columns.tsv` maps final columns back to raw FAMSA columns, both
1-based. Changing only trimAl mode or post-alignment QC thresholds reuses the
raw alignment; changing `min_protein_length` affects extraction as well as
trimming, and therefore also recomputes the affected alignments and trees.
BUSCO marker selection is independent of this sequence-length setting.
After sequence/alignment QC, the merge records coverage for each retained gene
tree using its actual tips and all selected species as the denominator. Coverage
is diagnostic: there is no overall or order-level coverage cutoff. Loci rejected
by other QC checks do not count toward species coverage.
The merge records each species' gene-tree count
in `species_coverage.tsv` and its distribution in `gene_trees.json`. The TSV's
`represented` column means only that the species appears in at least one retained
tree; it is not a statistical quality flag. ASTRAL refuses to run if any selected
species appears in zero retained trees, or no trees remain. The diagnostic files
are preserved on this failure, including `unrepresented_species` in the JSON.
There is no ten-gene stopping condition. Strongly uneven missing data can still
weaken relationships between clades even when these checks pass.

The former `phylogeny.min_occupancy` setting is removed from both marker selection
and gene-tree QC. Delete it from local configuration; an obsolete setting is
rejected rather than silently ignored. Ranking by coverage and limiting the
number of loci control selection without a fixed coverage cutoff. Removing the
cutoff does not guarantee sufficient overlap between clades. This is separate
from `selection.busco_threshold`, which selects species using their BUSCO summary
scores and is unchanged.

The obsolete settings `min_sites`, `min_informative_sites`, and
`min_markers_per_species` are rejected; delete them from local configuration.
The 100-site cutoff was redundant with the default 100-known-residue filter
per sequence. The ten-informative-site and ten-gene cutoffs had not been
validated for these data. Counts remain available for assessing information
content without imposing these numerical exclusion rules.

VeryFastTree 4.0.5 runs in double precision with `-lg -gamma`. Topology search
uses LG+CAT; `-gamma` rescales branch lengths and evaluates Gamma20 likelihoods.
It is not a full LG+Gamma topology search. Gene trees keep SH-like local support;
these values are not bootstrap proportions and no bootstrap cutoff is applied.
The species tree uses one estimated tree per locus, not bootstrap replicates.

ASTRAL-IV estimates topology, local posterior probabilities, and CASTLES-II
branch lengths from these gene trees. It receives the actual mean retained
alignment length, the root, and a fixed seed. Its output includes terminal
branch lengths in substitution/site units; this file is **not a time tree**.
Gene-tree estimation errors and model misspecification can affect both topology
and lengths. The workflow does not claim that approximate gene trees eliminate
these errors.

Alignments and gene trees default to four threads and an 8 GB scheduling budget
each. trimAl uses one thread and a 4 GB budget. ASTRAL defaults to 32 threads
and 64 GB. These are scheduling reservations,
not measured upper bounds or hard process memory limits. Snakemake limits the
sum of concurrent reservations; Slurm enforces the allocation's total memory.
Per-stage benchmarks report observed memory and time. Input preparation uses
one sequential pass per species during extraction and transposes the extracted
proteins once, avoiding one full proteome scan per marker. It never builds a
supermatrix. There is no automatic tool-comparison analysis.

## Optional time tree

Absolute dates require independently justified node-age calibrations. Neither
BUSCO results nor an NCBI taxonomy constraint tree provides ages. Provide a TSV
with these columns (the numbers below are **synthetic format examples**):

```text
taxa	min_age_ma	max_age_ma	source
Species_A,Species_B	90	110	Citation or calibration record
Species_C,Species_D	30	40	Citation or calibration record
```

`taxa` defines the MRCA of at least two exact tip labels in the inferred tree.
Bounds are positive ages in millions of years before present; equal bounds
specify an exact date. Sources are mandatory provenance. Check that each MRCA
corresponds to the intended biological node. Overlapping definitions of the same
node, contradictory ancestral bounds, and unknown species are rejected.

```yaml
phylogeny:
  dating:
    enabled: true
    calibrations: config/calibrations.tsv
```

Run the same launcher with target `timetree`; enabling dating also adds it to
`phylogeny`. There is no default fossil age and no implicit conversion from
root height 1 to millions of years. All selected species are treated as extant
with age zero. The wrapper converts positive Ma ages to LSD2's dates increasing
toward the present: a 90–110 Ma bound becomes `b(-110,-90)`, with every tip at 0.
MRCA expressions identify calibrated nodes, including the root.

Dating uses **LSD2 least squares with one estimated substitution rate**, weighted
by input branch lengths by default. It fits node dates on the existing rooted
CASTLES tree without a new alignment, concatenation, topology search, or root
search. The default estimates a common rate; lineage-rate variation can affect
its suitability and the resulting ages. It does not use treePL's penalized
likelihood, rate smoothing, cross-validation or stochastic optimization restarts.

```yaml
phylogeny:
  dating:
    mem_gb: 4
    lsd2:
      variance: 1
      variance_parameter: null
      numsites: null
```

`variance: 0` selects unweighted least squares, `1` weights using input lengths,
and `2` reweights after an initial fit. `variance_parameter` optionally sets
LSD2's positive variance offset (`-b`); null leaves its native automatic choice.
The standalone tool is serial; the rule fixes its thread count to 1 internally.
The 4 GB memory value
is a scheduling reservation, not a measured upper bound for the real dataset.
The dating stage does not use `phylogeny.seed` and does not simulate confidence
intervals. Rate partitions are not implemented by this wrapper.

`numsites` defaults to `total_gene_sites` in species-tree provenance: the sum of
trimmed alignment lengths for the genes actually retained for ASTRAL. With the
automatic variance offset, LSD2 uses the larger of the median branch length and
`10 / numsites`. This site count does **not** establish the effective sample size
of CASTLES lengths with missing species and variation in locus rates. An explicit
positive integer supports sensitivity analysis; it is recorded and never inferred
from the mean gene length. No concatenated alignment is constructed. An older
tree provenance file lacking `total_gene_sites` requires regeneration or a
justified explicit override.

The native invocation uses `-l -1 -u 0 -U 0` to retain zero-length internal
branches and permit zero time branches. It leaves the input substitution tree
unchanged. Each completed result must have a finite rate and objective, preserve
the rooted topology and species, and satisfy positive calibration bounds with
all tips at age zero. Missing, nonfinite, incomplete or inconsistent results fail
the stage; no fallback dating engine is used.

LSD2 prints native dates and lengths with six significant digits. The wrapper
reads node dates from the native time Nexus, verifies that its branches agree
within that rounding precision, and intersects the date rounding intervals with
hard calibrations and ancestor-age bounds. It fails if these cannot be satisfied.
It then writes parent-minus-child ages as branches with 17 significant digits,
so independent branch rounding does not accumulate into unequal root-to-tip
lengths. Every adjusted node age is recorded. The final Newick is independently
checked for ultrametricity, topology and calibration agreement.

With bounded calibrations alone, LSD2 can return a range of equally optimal
rates and dates instead of identifying a unique absolute scale. In that case,
the wrapper retains the native midpoint date tree, records the feasible intervals
and native warning, and marks `review_status: needs_review`. These feasible
intervals, including native annotations named `CI_date`, are **not statistical
confidence intervals**. This workflow does not request LSD2 confidence simulations.
Any native warning or rate at the tool's lower bound also triggers review.

The published `species_tree.dated.nwk` is in Ma. The retained
`lsd2.dated.date.nexus` is native time output, whereas `lsd2.fitted.nwk` contains
fitted **substitution lengths**. Stable internal node names are restored by
rooted-clade matching. `lsd2_runs/` retains input files, native results, commands
and logs even after a failed job. Published outputs are written only after
validation, with `provenance.json` last as the completion record.

Remove `dating.treepl`, `dating.command`, and `dating.threads` from old overrides.
The dating executable and its thread count are fixed by the rule. Obsolete settings fail
explicitly. Old CV/restart files from an earlier analysis are historical outputs;
the current rule's outputs and provenance identify the LSD2 run. Successful
inference stages remain reusable when only dating settings or outputs change.

Review calibration placement and sensitivity to defensible alternative bounds,
variance settings and site counts before scientific interpretation. Calibration,
gene-tree and branch-length uncertainty is not propagated into confidence
intervals. [The evaluation](dating_evaluation.md) records synthetic runtime
measurements and their limits.

## Retrieving TimeTree secondary calibrations with nwkit

```yaml
phylogeny:
  dating:
    calibration_source: timetree
    timetree:
      max_representatives: 64
      max_queries: 32
      min_studies: 5
      offline: false
```

```bash
# Infer/reuse the species tree and retrieve calibration candidates, without dating.
./run_pipeline.sh --software-deployment-method conda \
  --configfile config/mydata.yaml config/phylogeny.local.yaml \
  --cores 32 --resources mem_gb=128 -- phylogeny_calibrations

# Infer/reuse the tree, retrieve/reuse the calibrations, then run LSD2.
./run_pipeline.sh --software-deployment-method conda \
  --configfile config/mydata.yaml config/phylogeny.local.yaml \
  --cores 32 --resources mem_gb=128 -- timetree
```

The target can run end to end. `dating.enabled: true` additionally includes
dating in `phylogeny`; it is not needed for the explicit `timetree` target.
`calibration_source: file` retains the original manually supplied TSV route.
The file route can use an edited copy of the candidate TSV after literature
review; do not edit generated workflow outputs in place.

Retrieval uses the separate [timetree environment](../workflow/envs/timetree.yaml),
with Bioconda nwkit 0.27.0 and a verified client-module hash.
It leaves the inference environment unchanged.
The rule uses `python` from this environment; the interpreter is not a
configuration option.

The adapter calls `nwkit.mcmctree._fetch_timetree_url`, the same HTTP client used
by `nwkit mcmctree --timetree ci`. It requests the documented **JSON endpoint**
and interprets the response locally; it does not run the entire nwkit automatic
Newick-annotation command. This preserves returned MRCA IDs, found/missing
taxids, and bibliographic records. The private API boundary checks the pinned
module hash and is tested against the real installed client; a new nwkit
implementation requires an explicit adapter review. It does not modify nwkit.

The steps are:

1. Resolve metadata taxids against the existing local NCBI SQLite snapshot.
   Resolve merged IDs and subspecies to species; omit unresolved mappings and
   aliases assigned to the same species ID. Write `taxa.tsv`. No fuzzy name
   matching or explicit higher-rank TimeTree retry is used.
2. Partition the rooted inferred tree by repeatedly splitting its largest
   clade, up to 64 representative slots. Pick the tip with the highest retained
   gene-tree count in each part, breaking ties by name. Missing taxids can leave
   slots empty. A user file `timetree.representatives` can instead specify exact
   tree labels, one per line. Both root children must be represented. A
   high-coverage tip is not necessarily well represented in TimeTree.
3. Consider internal nodes with representatives on every child lineage, without
   a minimum number of original-tree tips. Query at most the 32 largest eligible
   clades, breaking size ties lexicographically by representative species labels.
   Each query uses the representative species IDs. The representative and query
   caps bound workload, not statistical accuracy. The former
   `phylogeny.dating.timetree.min_clade_taxa` setting is removed and rejected if
   supplied. Small clades, including two-tip clades, are eligible, but the ranking
   still favors larger clades and does not ensure even calibration placement.
4. Fetch serially with a one-second delay between uncached queries and nwkit's
   timeout/retry handling. Cache raw responses, URLs, timestamps, checksums and
   nwkit provenance. Network failures stop the stage; they are not converted
   into missing biological evidence. A successful API report of no estimate is
   retained as a missing result.
5. Require observed species on every input child lineage and check that their
   MRCA is still the intended input node. Require finite positive age ranges
   containing the point estimate, with at least five distinct bibliographic
   records. Reject point-only estimates and zero bounds. All assignments of a
   repeated TimeTree MRCA ID are excluded as ambiguous; equal ages alone are
   not treated as a duplicate. This conservative choice can lose usable dates.
6. Check ancestral/descendant bounds jointly. If they conflict, or no bounds
   survive, preserve diagnostics and write an empty calibration TSV so dating
   cannot proceed with an invalid set. Inspect `provenance.json` for `ready`,
   `conflicting_calibrations`, or `no_valid_calibrations`.

`candidates.tsv` and `candidates.json` describe each attempted node, its original
clade size, queried and used representatives, MRCA ID, study count, bounds and
exclusion reason. Raw response records in `resources/timetree_cache/` include `study_data`.
`representatives.nwk` is a pruned diagnostic copy; the original CASTLES tree
is the input to LSD2 and is never replaced by a calibration-only tree.

The cache has no automatic expiry: repeat runs use the frozen responses.
`offline: true` prohibits requests and fails on a cache miss. For an intentional
refresh, archive `resources/timetree_cache/` with the previous analysis record,
then rerun online with a new `analysis` name or explicitly force
`prepare_timetree_calibrations` using `--forcerun`. Removing the cache alone does
not invalidate completed calibration outputs. The fixed cache directory is
created as needed. Do not publish retrieved TimeTree datasets in this source
repository.

These are **secondary calibrations**, not fossil minima/maxima. TimeTree's
reported intervals summarize variation across published estimates, and fewer
studies can instead produce min–max ranges. They are not independent fossil
priors. This pipeline explicitly interprets retained bounds as hard LSD2
constraints; it does not transfer MCMCtree's soft tails or propagate a 95%
credible interval to the dated tree. Five papers need not be five independent
datasets. The sampled-MRCA checks do not prove that the full TimeTree topology
matches the inferred tree or that the intended crown/stem interpretation is
correct. Major calibrations need biological and original-study review.
See the [TimeTree FAQ](https://timetree.org/faqs),
[API examples](https://timetree.org/api), and
[nwkit calibration documentation](https://github.com/kfuku52/nwkit/wiki/nwkit-mcmctree).

## Outputs and validation

Outputs are under `results/<analysis>/phylogeny/` for `all`, and
`results/<analysis>/phylogeny_phenotyped/` for `phenotyped`:

| Output | Meaning |
| --- | --- |
| `selection/samples.tsv`, `selection/selection.json` (phenotyped only) | Observed-species input rows, trait column/source checksum, species counts and missing-trait exclusions |
| `plan/marker_stats.tsv`, `plan/markers.tsv` | All BUSCO marker statistics/ranks and the selected set in coverage-descending, ID-ascending order |
| `plan/species.tsv`, `plan/provenance.json` | Source paths, table checksums, lineage and selection settings |
| `species/*.faa`, `species/*.json` | Proteins; cdskit source hashes, padding/frame changes, masked codon positions and sequence QC |
| `markers/`, `alignments/raw/` | Per-marker inputs and reusable raw FAMSA alignments/QC |
| `alignments/*.faa`, `*.json`, `*.columns.tsv` | trimAl/QC results, site/variable/informative counts, mode, tool hash and final-to-raw column correspondence |
| `gene_trees/` | Per-marker trees and QC |
| `species_coverage.tsv`, `gene_trees.json` | Retained loci per species, presence flags, count distribution, absent species and exclusions |
| `species_tree.nwk`, `species_tree.json` | Rooted ASTRAL/CASTLES-II tree in substitutions/site |
| `taxonomy_audit/` | Optional MonoPhy species-tree review: groups, intruder/outlier flags, run linkage, rank figures and provenance |
| `timetree/calibrations.tsv`, `timetree/candidates.tsv`, `timetree/provenance.json` | Optional retrieved secondary bounds and mapping/exclusion audit |
| `dating/species_tree.dated.nwk`, `dating/node_ages.tsv` | Ultrametric time tree and node ages in Ma, with native ages and rounding adjustments |
| `dating/calibrations.resolved.tsv`, `dating/rounding_adjustments.tsv` | Resolved MRCA bounds and every age adjustment within native rounding precision |
| `dating/lsd2.dated.date.nexus`, `dating/lsd2.fitted.nwk` | Native time Nexus and fitted substitution tree, respectively |
| `dating/lsd2.report.txt`, `dating/lsd2.dates.txt`, `dating/lsd2.input.nwk`, `dating/lsd2.command.json` | Native fit report, date constraints, input tree and invocation |
| `dating/provenance.json`, `dating/lsd2_runs/` | Settings, units, fit/rate diagnostics, executable/build hashes and retained native runs/logs |
| `benchmarks/` | Runtime and maximum resident memory measurements |

Unit tests cover coverage-based selection and deterministic ties independently
of order labels, row order, replicate counts and BUSCO match lengths, retention
of loci below 50% coverage at selection and after QC, coverage diagnostics,
BUSCO lineage/ID handling, duplicate rejection, reading frame,
stop codons, ambiguous/partial-gap codons, alternate genetic codes, alignment
filtering, variable loci with low/zero informative-site counts, species with
only one retained gene tree, missing-species rejection, and calibration consistency.
cdskit API results are checked against
its actual CLI. Real-tool tests exercise FAMSA, trimAl, VeryFastTree, ASTRAL-IV and
LSD2, including column correspondence, original-X retention, unit conversion,
hard calibrations, native rounding and zero branches under all three variance
settings. They check unchanged reruns, dating-only recovery after a missing output,
and reuse of FAMSA after a trimAl mode change. Install the pinned
cdskit in the test Python environment. Enable the other tools with `SNAKEMAKE_BIN`,
`FAMSA_BIN`, `TRIMAL_BIN`, `VERYFASTTREE_BIN`, `ASTRAL_BIN`, and `LSD2_BIN` when they are
not on `PATH`. These variables are for tests, not production configuration.
The workflow integration requires the official ASTRAL build at its fixed project
location; the test exposes the other tools under their standard command names
on `PATH`. Setting `PHYLOGENY_CONDA_PREFIX` also runs the workflow integration
with actual stage-specific Conda environments. Synthetic tests establish workflow
behavior, not biological accuracy or 6,000-species performance.

`tests/test_phenotyped_phylogeny.py` checks observed zeros, missing annotations,
single-state and continuous traits, species counts across replicate rows, and
invalid set choices. Its real-tool workflow test runs phenotyped species first,
then all species, then both; verifies unchanged inference files across switches;
checks independent roots, calibrations and dating; and verifies that phenotype
membership changes leave the all-species outputs untouched. TimeTree responses
are recorded synthetic fixtures; these tests make no external API requests.

TimeTree tests also retain valid calibrations from small clades (including
two-tip clades) while checking that the query cap and size ranking still apply.
These tests use recorded synthetic responses without contacting TimeTree.

A separate real-input smoke check used Amborella trichopoda, Oryza sativa,
Abelia chinensis and Abeliophyllum distichum from tlight, with 20 shared Complete
BUSCO markers (sorted by ID; each match at least 100 aa). All 80 CDS passed
preparation, with no reading-frame changes; terminal stops were masked. FAMSA
and the previous trimAl/QC retained 18 markers and excluded two with eight
parsimony-informative sites each. With the simplified QC, all 20 are retained:
the two recovered alignments have 189/150 sites, including 134/113 variable
sites, respectively. The original 18 retained alignments and column maps are
unchanged. Original results are in `results/phylogeny_cdskit_trimal_smoke/summary.json`;
the repeated trimAl/QC and counts are in
`results/phylogeny_qc_simplification_smoke/summary.json`. The smoke check's
preselection by match length is not a rule in the production marker selection.
This checked input
compatibility through trimming; it did not infer a species tree or assess
full-scale accuracy/performance.

## Methods and source documentation

- [GeneGalleon species-tree workflow](https://github.com/kfuku52/genegalleon/blob/main/workflow/core/gg_genome_evolution_core.sh)
- [BUSCO output formats](https://busco.ezlab.org/busco_userguide.html)
- [FAMSA](https://github.com/refresh-bio/FAMSA)
- [cdskit preparation source](https://github.com/kfuku52/cdskit/tree/0.27.0/cdskit)
- [trimAl v1.5.1 source](https://github.com/inab/trimal/tree/d637091abe33595775f40480970d1a18d87a7bcb/source)
- [VeryFastTree](https://github.com/citiususc/veryfasttree)
- [FastTree models, support and Gamma20 scaling](https://morgannprice.github.io/fasttree/)
- [ASTRAL-IV and CASTLES-II, including the >5,000-species build requirement](https://github.com/chaoszhang/ASTER/blob/master/tutorial/astral4.md)
- [LSD2 source and dating options](https://github.com/tothuhien/lsd2)
- [LSD method paper](https://doi.org/10.1093/sysbio/syv068)
- [Coalescent branch lengths and scalable dating](https://doi.org/10.1093/sysbio/syag038)
