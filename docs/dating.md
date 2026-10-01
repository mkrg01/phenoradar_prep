# Dating species trees

[Documentation](index.md) · [BUSCO phylogeny](phylogeny.md)

LSD2 estimates ages in millions of years (Ma) using TimeTree calibrations or
manual bounds. It keeps the inferred topology and selected basal edge, while
estimating the root's position along that edge. All tips are treated as living
samples with age zero.

Before submitting an analysis, enable dating in `analysis.yaml`:

```yaml
phylogeny:
  trees: [all]
  dating:
    enabled: true
    calibration_source: timetree
```

Dating supports `all`, `phenotyped`, or both; representative trees are unsupported.

## TimeTree calibrations

Prepare calibrations for review; missing tree inference is included:

```bash
./run_analysis.sh submit --analysis results/leaf/downstream/carnivory --target phylogeny_calibrations
```

In `phylogeny/<set>/timetree/`, inspect `calibrations.tsv` for bounds,
`candidates.tsv` for selection decisions, and `studies.tsv` for sources.
`provenance.json` must report `status: ready`.

Candidates require every child lineage to be covered and at least five distinct
named study records. Accepted intervals become hard bounds, so review the
ancestors and studies; study count alone does not establish quality.
Missing or conflicting calibrations stop dating.

After review:

```bash
./run_analysis.sh submit --analysis results/leaf/downstream/carnivory --target timetree
```

`timetree` also prepares missing calibrations; enabled dating runs with `all`.
TimeTree requests deduplicate species taxids. Nodes whose child lineages share
a biological species, or have conflicting species-to-taxid assignments, are
excluded from calibration. To refresh cached responses, archive
`resources/timetree_cache/` and prepare a new analysis.

## Manual calibrations

Set `inputs.calibrations: input/calibrations.tsv` and
`phylogeny.dating.calibration_source: file`. Each row identifies an ancestor
by at least two exact sample tip labels, comma-separated:

```tsv
taxa	min_age_ma	max_age_ma	source
Species_A_run1,Species_B_run1	90	110	Citation or calibration record
```

The tips' most recent common ancestor receives the bounds. Use appropriate
samples, ages, and citations for your data. Bounds must be positive with
`min_age_ma <= max_age_ma`; equal bounds fix an age. Unknown tips, duplicate
ancestors, and conflicting bounds are rejected. Every row must apply to every
requested tree set.

## Check the results

Under `results/<build>/downstream/<analysis>/phylogeny/<set>/dating/`:

| File | Use |
| --- | --- |
| `species_tree.dated.nwk`, `node_ages.tsv` | Dated tree and ages in Ma |
| `provenance.json` | Settings, checks, solution bounds, and diagnostics |
| `calibrations.resolved.tsv` | Applied age bounds and node labels |
| `lsd2.input.nwk`, `lsd2.dates.txt`, `lsd2.outgroups.txt` | Native inputs |
| `lsd2.command.json`, `lsd2.log`, `lsd2.result` | Executed command and reports |
| `lsd2.result.date.nexus`, `lsd2.dated.nwk` | Native dated outputs |

`lsd2.result.nwk` is the fitted substitution tree, not the dated tree. Failed-job
logs remain in `lsd2_runs/`.

With living tips and interval calibrations, the absolute time scale can be
nonunique. LSD2 then exports the midpoint of its boundary solutions. Review
`unique_time_scale`, `root_age_solution_bounds_ma`, and `needs_review`.
**Solution bounds are not statistical confidence intervals**, even when native
NEXUS labels use `CI_*`. Do not replace bounds with exact ages just to remove
this diagnostic.

`checks_passed` confirms numerical validity, not biological accuracy. Ages depend
on the input tree, branch lengths, rate model, and calibrations; gene-tree and
calibration uncertainty are not propagated.

## Rooting and numerical checks

The [selected basal group](phylogeny.md#rooting) is retained. LSD2 optimizes only
that root edge (`-g`, `-r k`), retaining short/zero branches with minimum time
lengths of zero. It does not collapse support or remove outlier tips.

The wrapper checks tips, rooted topology, nonnegative lengths, ultrametricity,
and calibration bounds. Tiny negative lengths at floating-point zero are set
to zero and recorded in `numerical_zero_adjustments.tsv`; larger negatives fail.
Raw native files preserve original values.

## Optional overrides

Normally omit `phylogeny.dating.lsd2`. Supported settings are:

| Setting | Default | Meaning |
| --- | --- | --- |
| `variance` | `1` | `0`: unweighted; `1`: input-length weights; `2`: refit with fitted-length weights |
| `variance_parameter` | `null` | Automatic variance offset, or a positive value up to 1 |

No branch-rate partitions or confidence-interval simulations are requested.
The old `phylogeny.dating.treepl` setting must be removed for new LSD2 analyses.

## Resources

`date_busco_species_tree` requires one CPU and defaults to 4 GB per tree set.
Adjust memory/time through [resource overrides](running.md#resource-budgets).

## Methods

[LSD2 source and options](https://github.com/tothuhien/lsd2) ·
[Least-squares dating](https://doi.org/10.1093/sysbio/syv068) ·
[Coalescent branch lengths and dating](https://doi.org/10.1093/sysbio/syag038)
