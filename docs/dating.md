# Dating species trees

[Documentation](index.md) · [BUSCO phylogeny](phylogeny.md)

The `timetree` target dates rooted BUSCO trees with LSD2 using manual or TimeTree
calibrations. Input branch lengths are substitutions/site; output ages are in
millions of years (Ma). The target follows `phylogeny.species_sets`.

## Manual calibrations

Provide independently justified node-age bounds in a TSV. These numbers are
synthetic format examples:

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

```bash
./run_pipeline.sh --configfile config/mydata.yaml \
  --cores 32 --resources mem_gb=128 -- timetree
```

Missing tree inference is scheduled automatically. The explicit target works
without `dating.enabled`; enabling it adds dating to `phylogeny`. All tips are
treated as extant at age zero. Calibrations must fit every requested species set.

## LSD2 fitting

LSD2 preserves the input topology and root, fitting a single rate with no rate
partitions or confidence simulations. The container supplies LSD2 2.4.4; native
setup uses [dating.yaml](../workflow/envs/dating.yaml) and its verified source
installer. `LSD2_SOURCE_ARCHIVE` can supply the pinned source offline.

```yaml
phylogeny:
  dating:
    mem_gb: 4
    lsd2:
      variance: 1
      variance_parameter: null
      numsites: null
```

`variance` selects unweighted fitting (`0`), input-length weights (`1`), or
reweighting after an initial fit (`2`). `variance_parameter` sets a positive
variance offset; null uses LSD2's automatic choice. The rule uses one thread.

`numsites` defaults to the sum of retained gene-alignment lengths recorded in
the species-tree provenance. It affects variance weighting and is not an
estimate of effective sample size. Older provenance lacking `total_gene_sites`
requires regeneration or an explicit positive override.

The wrapper validates finite results, topology, bounds, and ultrametricity,
and reconciles native rounding before writing the final time tree. Native files
and failed runs remain in `dating/lsd2_runs/`.

Review `dating/provenance.json`, especially `review_status: needs_review` and
native warnings. Bounded calibrations can leave a range of equally optimal ages;
reported feasible ranges, including native `CI_date`, are not statistical
confidence intervals. Assess sensitivity to calibrations and variance settings.
Gene-tree and calibration uncertainty are not propagated into confidence intervals.

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
./run_pipeline.sh --configfile config/mydata.yaml \
  --cores 32 --resources mem_gb=128 -- phylogeny_calibrations
```

This infers/reuses trees and retrieves candidates without dating. Review
`timetree/calibrations.tsv` and `candidates.tsv`, then run `timetree`. To edit
bounds manually, copy the TSV to `input/` and select `calibration_source: file`.

Retrieval resolves NCBI species taxids, partitions the tree into representative
clades, and favors tips with more retained gene trees. An optional
`timetree.representatives` file lists exact tree labels, one per line; both root
children must be covered. Eligible query nodes need representatives on every
child lineage; larger clades have priority within the query cap.

Accepted estimates need valid positive bounds, the intended sampled MRCA, and
at least `min_studies` distinct records. Repeated TimeTree MRCA assignments are
excluded. Conflicting or absent bounds leave an empty calibration TSV and a
readiness explanation in `timetree/provenance.json`; dating cannot proceed.

Responses are cached in `resources/timetree_cache/` without automatic expiry.
`offline: true` forbids requests and fails on a cache miss. Network errors stop
retrieval. To refresh deliberately, archive the cache and rerun online with a
new `run_name` or `--forcerun prepare_timetree_calibrations`.

TimeTree estimates are secondary calibrations summarizing published studies.
They become hard LSD2 bounds, not fossil priors or confidence intervals for the
resulting tree. Review key MRCA placements and original studies; see the
[TimeTree FAQ](https://timetree.org/faqs) and
[nwkit documentation](https://github.com/kfuku52/nwkit/wiki/nwkit-mcmctree).

## Outputs

Under `results/<run_name>/phylogeny/<set>/`:

| Output | Contents |
| --- | --- |
| `timetree/calibrations.tsv`, `candidates.tsv`, `provenance.json` | Retrieved bounds, rejected candidates, readiness status |
| `timetree/taxa.tsv`, `representatives.nwk` | Taxid resolution and representative subtree |
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
