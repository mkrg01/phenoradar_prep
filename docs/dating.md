# Dating species trees

[Documentation](index.md) · [BUSCO phylogeny](phylogeny.md)

The `timetree` target retrieves TimeTree calibrations by default and dates rooted
BUSCO trees with LSD2. Input branch lengths are substitutions/site; output ages are in
millions of years (Ma). The target follows `phylogeny.species_sets`.

```bash
./run_pipeline.sh --configfile config/mydata.yaml \
  --cores 32 --resources mem_gb=128 -- timetree
```

Missing tree inference is scheduled automatically. The explicit target works
without `dating.enabled`; enabling it adds dating to `phylogeny`. All tips are
treated as extant at age zero. Calibrations must fit every requested species set.

## Manual calibrations

To use your own bounds, set `calibration_source: file` and provide a TSV.
These numbers are synthetic format examples:

```text
taxa	min_age_ma	max_age_ma	source
Species_A,Species_B	90	110	Citation or calibration record
Species_C,Species_D	30	40	Citation or calibration record
```

`taxa` identifies the MRCA of at least two exact tree labels. Bounds must be
positive Ma; equal bounds specify an exact age. `source` is required. Check the
biological placement of each MRCA; unknown tips and conflicting bounds are rejected.

```yaml
phylogeny:
  dating:
    enabled: true
    calibration_source: file
    calibrations: input/calibrations.tsv
```

## LSD2 fitting

LSD2 preserves the input topology and root, fitting a single rate with no rate
partitions or confidence simulations. The container supplies LSD2 2.4.4; native
setup uses [dating.yaml](../workflow/envs/dating.yaml) and its verified source
installer. `LSD2_SOURCE_ARCHIVE` can supply the pinned source offline.

Fitting always uses input-length weights (`-v 1`) and LSD2's automatic variance
offset (no `-b` override). The site count (`-s`) is the sum of retained trimmed
gene-alignment lengths from `total_gene_sites` in the species-tree provenance.
It affects variance weighting and is not an estimate of effective sample size.
Older provenance lacking this count requires regeneration. These choices are
fixed in the dating helper and require no configuration.

The `date_busco_species_tree` rule defaults to 1 thread and 4 GB;
see [resource overrides](running.md#resource-budgets).

The wrapper validates finite results, topology, bounds, and ultrametricity,
and reconciles native rounding before writing the final time tree. Native files
and failed runs remain in `dating/lsd2_runs/`.

Review `dating/provenance.json`, especially `review_status: needs_review` and
native warnings. Bounded calibrations can leave a range of equally optimal ages;
reported feasible ranges, including native `CI_date`, are not statistical
confidence intervals. Assess sensitivity to calibrations.
Gene-tree and calibration uncertainty are not propagated into confidence intervals.

## TimeTree calibrations

```yaml
phylogeny:
  dating:
    calibration_source: timetree
```

```bash
./run_pipeline.sh --configfile config/mydata.yaml \
  --cores 32 --resources mem_gb=128 -- phylogeny_calibrations
```

Use this optional step to inspect retrieved candidates before dating. Review
`timetree/calibrations.tsv`, `candidates.tsv`, and `studies.tsv`, then run `timetree`. To edit
bounds manually, copy the TSV to `input/` and select `calibration_source: file`.

Retrieval follows a fixed procedure with no TimeTree configuration section:

1. Resolve the full tree's species to NCBI species taxids using the local taxonomy
   snapshot. Ambiguous mappings and labels sharing a species taxid are excluded
   from queries; those species remain in the inference and dating tree.
2. Consider every internal node. Query all resolved descendant species when
   every child lineage has at least one. There is no representative selection,
   query-count cap, or clade-size threshold. Visit nodes in canonical preorder,
   ordering child clades by their first species label. Query taxids are sorted.
3. Require returned taxa on every child lineage, a returned MRCA identifier, and
   finite bounds satisfying `0 < lower <= reported age <= upper` with
   `lower < upper`. Require at least five distinct named bibliographic records
   returned for that MRCA. Count title/year/first-author combinations, ignoring
   case and repeated whitespace; multiple trees from a paper do not add studies.
   The five-study condition is fixed in the helper. Returned intervals are used
   as supplied, without recalculating them from the bibliographic list.
4. Exclude all ambiguous assignments when the same returned MRCA identifies
   multiple input nodes. Check that the remaining bounds permit ancestors to
   be at least as old as their descendants. Conflicting or absent calibrations
   leave an empty calibration TSV and a readiness explanation in
   `timetree/provenance.json`; dating stops. Bounds are never widened and a
   compatible subset is never chosen automatically.

Every internal node has a row in `candidates.tsv`, including nodes that could
not be queried. Rows record query/observed taxids, the returned MRCA, study count,
available ages, exclusion reason, URL, and original retrieval time. `studies.tsv`
contains the deduplicated bibliographic records for both retained and excluded
candidates, including each record's original JSON. Full API responses remain in
the cache, and `candidates.json` records their paths and checksums. `nodes.nwk`
retains all species and substitution lengths with stable internal labels linking
the reports to clades. The input species tree is unchanged.

Responses are cached in `resources/timetree_cache/` without automatic expiry.
Cached responses are reused; only missing responses are retrieved from TimeTree.
No network connection is needed when all required responses are cached. Network
errors stop retrieval. To refresh deliberately, archive the cache and rerun with
a new `run_name` or `--forcerun prepare_timetree_calibrations`.

TimeTree estimates are secondary calibrations summarizing published studies.
They become hard LSD2 bounds, not fossil priors or confidence intervals for the
resulting tree. TimeTree describes its interval for five or more studies as
among-study variation under the 95% Empirical Rule; fewer studies receive a
min/max range. Our count uses the bibliographic records supplied by the API,
not an independently reconstructed set of observations behind that interval.
Five records do not establish statistical independence or calibration quality.
Child-lineage coverage and distinct MRCA identifiers do not establish full
topological concordance with TimeTree. Review key MRCA placements and original studies; see the
[TimeTree FAQ](https://timetree.org/faqs) and
[nwkit documentation](https://github.com/kfuku52/nwkit/wiki/nwkit-mcmctree).

### Reporting the method

For a completed run, a Methods description can follow this template. Fill in
the versions, retrieval dates, and retained calibration count from provenance:

> Divergence times were estimated on the rooted species-tree topology using
> substitution-per-site branch lengths and weighted least squares in LSD2
> [version]. All extant tips were assigned age zero. All internal nodes were
> considered for secondary calibration using TimeTree responses retrieved on
> [dates] through the nwkit [version] HTTP client. We retained positive age
> intervals associated with at least five distinct bibliographic records and
> observed taxa spanning every child lineage of the target node. Ambiguous
> repeated assignments of a TimeTree MRCA were excluded. The [N] retained
> intervals were checked for temporal compatibility and imposed as hard bounds.
> The topology and root were fixed, and a single substitution rate was fitted
> with input-length weights. Calibration records and exclusion reasons were
> documented in the supplementary tables.

The report is a point estimate conditional on these calibrations. For a paper,
assess sensitivity by omitting each calibration in turn where the remaining
calibrations still permit dating, then checking important ages and downstream
conclusions. This additional analysis is not run automatically.

## Outputs

Under `results/<run_name>/phylogeny/<set>/`:

| Output | Contents |
| --- | --- |
| `timetree/calibrations.tsv`, `candidates.tsv`, `provenance.json` | Accepted bounds, every candidate and its disposition, readiness status |
| `timetree/candidates.json`, `studies.tsv` | Cache checksums and bibliographic records linked to candidate nodes |
| `timetree/taxa.tsv`, `nodes.nwk` | Taxid resolution and full tree with report node labels |
| `dating/species_tree.dated.nwk`, `node_ages.tsv` | Time tree and node ages in Ma |
| `dating/calibrations.resolved.tsv`, `rounding_adjustments.tsv` | Resolved bounds and rounding adjustments |
| `dating/lsd2.dated.date.nexus` | Native time output |
| `dating/lsd2.fitted.nwk` | Fitted substitution lengths, not ages |
| `dating/provenance.json`, `lsd2_runs/` | Fit diagnostics, settings, and native runs/logs |

Dating-only changes reuse completed inference. Logs and benchmarks follow the
same branch under `logs/<run_name>/phylogeny/`.

## Method references

- [LSD2 source and options](https://github.com/tothuhien/lsd2)
- [LSD method paper](https://doi.org/10.1093/sysbio/syv068)
- [Coalescent branch lengths and scalable dating](https://doi.org/10.1093/sysbio/syag038)
