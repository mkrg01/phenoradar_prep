# Dating species trees

[Documentation](index.md) · [BUSCO phylogeny](phylogeny.md)

`timetree` uses **LSD2** to estimate ages in millions of years (Ma), with
**TimeTree** calibrations by default or manual age bounds. The ASTRAL-IV topology
is fixed. LSD2 estimates the root position along the edge separating the selected
basal group from the remaining samples. All basal species and samples are retained;
all tips are treated as living samples with age zero.

Before [preparing analysis](datasets.md#run-an-analysis), set
`phylogeny.trees: [all]`, `[phenotyped]`, or both, and
`phylogeny.dating.enabled: true`. Representative trees are unsupported.

## TimeTree calibrations

Prepare and review calibrations before dating; missing tree inference is included:

```bash
./run_analysis.sh submit --analysis analyses/analysis001 --target phylogeny_calibrations
```

Inspect `results/<analysis>/phylogeny/<set>/timetree/`: `calibrations.tsv` has bounds,
`candidates.tsv` explains selection, and `studies.tsv` lists sources. Check
`provenance.json` for `status: ready` and review the ancestors and supporting studies.

Candidates require coverage of every child lineage and at least five distinct
named study records. Accepted intervals become hard bounds; study count alone
does not establish quality. Missing/conflicting calibrations stop dating.

After that job finishes and the results are reviewed:

```bash
./run_analysis.sh submit --analysis analyses/analysis001 --target timetree
```

`timetree` prepares missing calibrations too. Enabled dating is included in `all`;
`phylogeny` stops at inference. TimeTree responses are cached across analyses.
To refresh, archive `resources/timetree_cache/` and prepare a new analysis.

Multiple samples sharing one biological `species_id` and taxid are retained.
TimeTree requests deduplicate taxids; nodes whose child lineages share a species
are excluded from calibration rather than treated as a species divergence.
Conflicting biological species mapped to the same species taxid remain excluded.

## Manual calibrations

Set `inputs.calibrations: input/calibrations.tsv` and
`phylogeny.dating.calibration_source: file` in `analysis.yaml` before preparing.
Each row identifies an ancestor by at least two exact sample tip labels,
comma-separated. Their MRCA after basal-group rooting receives the bounds:

```tsv
taxa	min_age_ma	max_age_ma	source
Species_A_run1,Species_B_run1	90	110	Citation or calibration record
```

Replace the example with appropriate samples, ages, and citations. Bounds must be
positive with `min_age_ma <= max_age_ma`; equal bounds fix an age. Unknown tips,
duplicate ancestors, and conflicting bounds are rejected. Calibrations must
apply to every requested species set. Ages in Ma are converted internally to
negative LSD2 dates: `[90,110] Ma` becomes `b(-110,-90)`, with tips at zero.

## Rooting and topology

Use `phylogeny.outgroup` to select a basal species or a list of basal species;
all their selected samples are included. The group must be separable by one edge
of the inferred molecular tree. Nonmonophyletic groups stop the analysis.
No species are added for rooting, and no basal tips are removed.

LSD2 receives the complete group via `-g` and optimizes only that edge (`-r k`).
It does not search all root edges. The two input root branches are treated as
one edge whose split is re-estimated; an input 50:50 split does not fix their
dated lengths. Short and zero-length branches are retained (`-l -1`), with
minimum time lengths zero (`-u 0 -U 0`). No support collapse or outlier-tip
removal is requested. The wrapper checks all tips, the rooted topology,
nonnegative lengths, ultrametricity and calibration bounds at native precision.

The source substitution tree remains available for contrast analysis. Its root
marks the basal split; its root-edge split is a serialization convention.
See [ASTRAL rooting](phylogeny.md#rooting) for its single-anchor limitation.
LSD2 does not correct upstream topology or substitution-length errors.

## Check the results

Under `results/<analysis>/phylogeny/<set>/dating/`:

| File | Use |
| --- | --- |
| `species_tree.dated.nwk`, `node_ages.tsv` | Dated tree and ages in Ma |
| `provenance.json` | Version, settings, root handling, solution bounds and diagnostics |
| `calibrations.resolved.tsv` | Applied bounds and pipeline node labels |
| `lsd2.input.nwk`, `lsd2.dates.txt`, `lsd2.outgroups.txt` | Reproducible native inputs |
| `lsd2.command.json`, `lsd2.log`, `lsd2.result` | Executed command and native reports |
| `lsd2.result.date.nexus`, `lsd2.dated.nwk` | Native time tree and extracted Newick |
| `lsd2.result.nwk` | Native fitted substitution tree; not the time tree |

Native logs also remain in `dating/lsd2_runs/` if a job fails. The pipeline tree
restores stable internal labels and native time-edge values, except for the
explicit numerical-zero handling below.
Tiny negative lengths at numerical zero (at most 64 floating-point ulps of the
root-age scale) are set to zero and listed in `numerical_zero_adjustments.tsv`;
larger negative lengths fail validation. Raw native files retain the original
values. Node ages are mean descendant-tip distances from the published edges.

With contemporaneous tips and only interval calibrations, the absolute time scale
can be nonunique. LSD2 then exports the midpoint of its boundary-date solutions.
`unique_time_scale: false`, `root_age_solution_bounds_ma`, and `needs_review`
identify this case. These solution bounds are **not statistical confidence
intervals**, even though native NEXUS labels use `CI_*`. Do not replace intervals
with exact ages solely to suppress this diagnostic.

`checks_passed` means numerical checks passed. Ages remain conditional on the
input tree, lengths, rate model and calibrations; gene-tree and calibration
uncertainty are not propagated. Native warnings appear in `diagnostics`.

## Optional overrides

Normally omit `phylogeny.dating.lsd2`. Supported options are:

| Setting | Default | Use |
| --- | --- | --- |
| `variance` | `1` | `0`: unweighted; `1`: weights from input lengths; `2`: second fit with weights from fitted lengths |
| `variance_parameter` | `null` | Native automatic variance offset; optionally supply a positive value up to 1 |

No branch-rate partitions or confidence-interval simulations are requested.
The retained site count is passed as `-s` for native variance defaults; it is
not treated as the effective sample size of a concatenated alignment.
Point dating is deterministic and takes no seed; upstream tree inference still
uses the analysis seed. Changed scientific settings require a new analysis.

The former `phylogeny.dating.treepl` configuration is no longer accepted.
Remove it when preparing an LSD2 analysis. Historical result directories are
left intact; new runs publish only LSD2 outputs as declared workflow products.

## Resources

`date_busco_species_tree` requires **one CPU** and defaults to **4 GB per species set**.
Check `logs/<analysis>/phylogeny/<set>/benchmarks/dating.tsv` and adjust memory/time
with [rule overrides](running.md#resource-budgets). Known basal-edge optimization
avoids an all-edge root search. Actual runtime depends on the tree and constraints.

## Method references

The dating environment builds unmodified LSD2 2.4.4 at commit
`c61110f3a4fa05325b45c97b2134792ff9d55d4c`; archive and executable checksums are
verified. Snakemake runs `workflow/envs/dating.post-deploy.sh` when deploying the
environment. For manual Conda installation, activate `dating.yaml`'s environment
and run that script; `LSD2_SOURCE_ARCHIVE` supports a verified offline archive.

[LSD2 source and options](https://github.com/tothuhien/lsd2) ·
[Least-squares dating method](https://doi.org/10.1093/sysbio/syv068) ·
[Coalescent branch lengths and dating](https://doi.org/10.1093/sysbio/syag038)
