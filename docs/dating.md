# Dating species trees

[Documentation](index.md) · [BUSCO phylogeny](phylogeny.md)

The `timetree` target dates the rooted BUSCO species tree with LSD2. Calibration
bounds come from a supplied TSV or retrieved TimeTree estimates. The output tree
uses millions of years (Ma); the input CASTLES tree uses substitutions per site.
Both `all` and `phenotyped` runs follow `phylogeny.species_sets`.

## Software

[`workflow/envs/dating.yaml`](../workflow/envs/dating.yaml) supplies Python and
compiler dependencies. Its [post-deploy script](../workflow/envs/dating.post-deploy.sh)
builds unmodified LSD2 2.4.4 from verified upstream commit
`c61110f3a4fa05325b45c97b2134792ff9d55d4c`, using Conda GNU C++ 13 and the upstream
makefile. Build provenance is in `share/lsd2/build.json` inside the environment.

For an offline source build, set `LSD2_SOURCE_ARCHIVE` to the verified archive
for that commit. Conda dependencies must also be available locally. TimeTree
retrieval uses the separate [nwkit environment](../workflow/envs/timetree.yaml).
The [tool evaluation](notes/dating_evaluation.md) records the installation choice
and synthetic benchmarks.

## Manual calibrations

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
    calibration_source: file
    calibrations: input/calibrations.tsv
```

Run the target `timetree`; enabling dating also adds it to `phylogeny`:

```bash
./run_pipeline.sh --software-deployment-method conda \
  --configfile config/mydata.yaml --cores 32 --resources mem_gb=128 -- timetree
```

This can schedule missing inference. Completed unchanged trees are reused.
Calibrations supply the absolute scale. All selected species are treated as extant
with age zero. The wrapper converts positive Ma ages to LSD2's dates increasing
toward the present: a 90–110 Ma bound becomes `b(-110,-90)`, with every tip at 0.
MRCA expressions identify calibrated nodes, including the root.

## LSD2 fitting

LSD2 fits node dates on the existing rooted CASTLES substitution tree, with
one estimated rate and input-length variance weights by default. It preserves
the topology and root. Lineage-rate variation can affect the suitability of
that common rate and the resulting ages.

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
The rule uses one thread and a 4 GB scheduling reservation. It does not use
`phylogeny.seed`. The wrapper supports one rate, without rate partitions or
confidence simulations.

`numsites` defaults to `total_gene_sites` in species-tree provenance: the sum of
trimmed alignment lengths for the genes actually retained for ASTRAL. With the
automatic variance offset, LSD2 uses the larger of the median branch length and
`10 / numsites`. This site count does **not** establish the effective sample size
of CASTLES lengths with missing species and variation in locus rates. An explicit
positive integer supports sensitivity analysis and is recorded. An older
tree provenance file lacking `total_gene_sites` requires regeneration or a
justified explicit override.

The native invocation uses `-l -1 -u 0 -U 0` to retain zero-length internal
branches and permit zero time branches. It leaves the input substitution tree
unchanged. Each completed result must have a finite rate and objective, preserve
the rooted topology and species, and satisfy positive calibration bounds with
all tips at age zero. Missing, nonfinite, incomplete, or inconsistent results fail the stage.

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

Review calibration placement and sensitivity to defensible alternative bounds,
variance settings and site counts before scientific interpretation. Calibration,
gene-tree and branch-length uncertainty is not propagated into confidence
intervals. [The evaluation](notes/dating_evaluation.md) records synthetic runtime
measurements and their limits.

## TimeTree calibrations

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
  --configfile config/mydata.yaml \
  --cores 32 --resources mem_gb=128 -- phylogeny_calibrations

# Infer/reuse the tree, retrieve/reuse the calibrations, then run LSD2.
./run_pipeline.sh --software-deployment-method conda \
  --configfile config/mydata.yaml \
  --cores 32 --resources mem_gb=128 -- timetree
```

`phylogeny_calibrations` requests TimeTree candidates. `timetree` uses
`calibration_source` to select retrieved or manual bounds, and works regardless
of `dating.enabled`. To review candidates manually, copy the generated TSV to
an input file, edit that copy, and select `calibration_source: file`.

Retrieval uses nwkit 0.27.0's verified HTTP client to fetch TimeTree JSON,
retaining MRCA IDs, taxids, and bibliography. Candidate selection proceeds as
follows (counts below are defaults):

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
   caps bound workload. Small clades, including two-tip clades, are eligible;
   larger clades retain priority, so calibration placement may be uneven.
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
not invalidate completed calibration outputs. Keep retrieved datasets with the analysis records, outside this source repository.

These are **secondary calibrations** summarizing variation among published
estimates; with fewer studies, reported bounds can be min–max ranges. The
workflow treats them as hard LSD2 constraints. They are not fossil priors or
confidence intervals for the resulting dated tree, and different papers may
share underlying data. Sampled-MRCA checks do not establish correspondence of
the full topologies or crown/stem interpretations. Review major calibrations
against the biological question and original studies.
See the [TimeTree FAQ](https://timetree.org/faqs),
[API examples](https://timetree.org/api), and
[nwkit calibration documentation](https://github.com/kfuku52/nwkit/wiki/nwkit-mcmctree).

## Outputs

Files are written under `results/<analysis>/phylogeny/<set>/`, where `<set>` is
`all` or `phenotyped`. Calibration labels must be valid for every requested tree;
automatic representatives and bounds are selected separately for each.

| Output | Contents |
| --- | --- |
| `timetree/taxa.tsv`, `representatives.nwk` | Taxid resolution and representative subtree |
| `timetree/calibrations.tsv`, `candidates.tsv`, `candidates.json`, `provenance.json` | Retrieved bounds, attempted queries, exclusion reasons, readiness status |
| `dating/species_tree.dated.nwk`, `node_ages.tsv` | Time tree and node ages in Ma |
| `dating/calibrations.resolved.tsv`, `rounding_adjustments.tsv` | Resolved MRCA bounds and native-rounding adjustments |
| `dating/lsd2.dated.date.nexus` | Native time output |
| `dating/lsd2.fitted.nwk` | Fitted substitution lengths; not a time tree |
| `dating/lsd2.report.txt`, `lsd2.dates.txt`, `lsd2.input.nwk`, `lsd2.command.json` | Native fit report, constraints, input tree, and invocation |
| `dating/provenance.json`, `lsd2_runs/` | Fit/rate diagnostics, units, settings, executable/build hashes, retained native runs and logs |

Logs and benchmarks use the matching branch under `logs/<analysis>/phylogeny/`.
Changes confined to dating settings or missing dated outputs reuse completed
inference. Review `dating/provenance.json` for warnings before using the tree.

## Method references

- [LSD2 source and options](https://github.com/tothuhien/lsd2)
- [LSD method paper](https://doi.org/10.1093/sysbio/syv068)
- [Coalescent branch lengths and scalable dating](https://doi.org/10.1093/sysbio/syag038)
