# Dating species trees

[Documentation](index.md) · [BUSCO phylogeny](phylogeny.md)

The `timetree` target dates rooted BUSCO trees with LSD2, using TimeTree
calibrations by default. Input lengths are substitutions/site; output ages are
millions of years (Ma). All tips are treated as extant at age zero.

```bash
./run_pipeline.sh --cores 32 --resources mem_gb=128 -- timetree
```

The target follows `phylogeny.species_sets` and schedules missing inference.
`phylogeny.dating.enabled: true` also adds dating to `phylogeny`.
Calibrations must fit every requested set.

## Manual calibrations

Provide a TSV with at least two exact tree labels per `taxa` entry, identifying
their most recent common ancestor (MRCA). Bounds must be positive Ma; equal
bounds specify an exact age. `source` is required. These are synthetic examples:

```text
taxa	min_age_ma	max_age_ma	source
Species_A,Species_B	90	110	Citation or calibration record
Species_C,Species_D	30	40	Citation or calibration record
```

```yaml
phylogeny:
  dating:
    calibration_source: file
    calibrations: input/calibrations.tsv
```

Check each MRCA's biological placement. Unknown tips and conflicting bounds
are rejected.

## TimeTree calibrations

With `phylogeny.dating.calibration_source: timetree` (default), inspect
calibrations before dating:

```bash
./run_pipeline.sh --cores 32 --resources mem_gb=128 -- phylogeny_calibrations
```

Review `timetree/calibrations.tsv`, `candidates.tsv`, and `studies.tsv` under each
species set, then run `timetree`. To edit bounds, copy the calibration TSV to
`input/` and select `calibration_source: file`.

Retrieval considers every internal node, querying its resolved descendant NCBI
species taxids. Ambiguous mappings and shared taxids are excluded from queries;
all species remain in the tree. A calibration is retained only when:

- Returned taxa span every child lineage of the target node.
- The returned MRCA has finite ages with `0 < lower <= age <= upper` and
  `lower < upper`, plus at least five distinct named bibliographic records.
  Multiple trees from one paper count once.
- Its returned MRCA does not ambiguously identify multiple input nodes.

Bounds must permit ancestors to be at least as old as descendants. Conflicting
or absent calibrations stop dating; `timetree/provenance.json` explains why.
The workflow never widens bounds or chooses a compatible subset automatically.

Responses are cached in `resources/timetree_cache/` without expiry. Only missing
responses need network access; network errors stop retrieval. To refresh, archive
the cache and rerun with a new `run_name` or
`--forcerun prepare_timetree_calibrations`.

TimeTree's secondary calibration intervals become hard LSD2 bounds as supplied.
Five records do not establish independence or quality; child-lineage coverage
does not ensure topological agreement. Review key MRCA placements and original
studies; see the [TimeTree FAQ](https://timetree.org/faqs).

## LSD2 fitting

LSD2 fixes topology/root and fits one rate with input-length weights (`-v 1`)
and an automatic variance offset. The site count is `total_gene_sites` from
species-tree provenance (summed retained alignment lengths), used for weighting,
not effective sample size. Older results lacking it need regeneration.

Review `dating/provenance.json` for `review_status: needs_review` and native
warnings. Validation checks topology, bounds, finite ages, and ultrametricity;
native rounding adjustments are recorded.

The tree is a point estimate conditional on calibrations. Feasible ranges,
including native `CI_date`, are **not statistical confidence intervals**;
gene-tree/calibration uncertainty is not propagated. Assess sensitivity by
omitting individual calibrations where dating remains possible; this is manual.

## Outputs

Under `results/<run_name>/phylogeny/<set>/`:

| Output | Use |
| --- | --- |
| `timetree/calibrations.tsv`, `provenance.json` | Accepted bounds and readiness status |
| `timetree/candidates.tsv`, `studies.tsv` | Every node's acceptance/exclusion reason and supporting studies |
| `timetree/taxa.tsv`, `nodes.nwk`, `candidates.json` | Taxid resolution, labeled clades, and cached-response checksums |
| `dating/species_tree.dated.nwk`, `node_ages.tsv` | Time tree and node ages in Ma |
| `dating/calibrations.resolved.tsv`, `rounding_adjustments.tsv` | Applied bounds and numerical adjustments |
| `dating/provenance.json`, `lsd2_runs/` | Fit diagnostics, settings, and native files/logs |

`dating/lsd2.fitted.nwk` contains fitted substitution lengths, not ages.
Dating-only changes reuse inference. For reporting, record tool versions,
retrieval dates, retained calibrations, fitting settings, and sensitivity results.

## Method references

- [LSD2 source and options](https://github.com/tothuhien/lsd2)
- [LSD method paper](https://doi.org/10.1093/sysbio/syv068)
- [Coalescent branch lengths and scalable dating](https://doi.org/10.1093/sysbio/syag038)
