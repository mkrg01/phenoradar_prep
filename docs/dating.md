# Dating species trees

[Documentation](index.md) · [BUSCO phylogeny](phylogeny.md)

The `timetree` target uses **treePL** to estimate divergence ages in millions of
years (Ma). **TimeTree** supplies age calibrations by default; you can also supply
your own minimum/maximum ages for ancestors. The topology and root stay fixed,
and all tips are treated as living species with age zero.

Prepare the [BUSCO inputs](phylogeny.md#inputs) and [software](containers.md), then
run from the repository root. `phylogeny.species_sets` selects the trees to date
(`[all]` by default). Missing upstream steps are run automatically.

## TimeTree calibrations

For a first run, prepare calibrations before dating:

```bash
./run_pipeline.sh --cores 32 --resources mem_gb=128 -- phylogeny_calibrations
```

Inspect `results/<run_name>/phylogeny/<set>/timetree/`: `calibrations.tsv` contains
age bounds, `candidates.tsv` explains selection, and `studies.tsv` lists sources.
Check that `provenance.json` reports `status: ready`, and review whether the
calibrated ancestors and supporting studies are appropriate for your analysis.

Candidates require coverage of every child lineage and at least five distinct
named study records. Missing or conflicting calibrations stop dating. Accepted
intervals become hard age bounds; the study count alone does not establish quality.

Run dating after review:

```bash
./run_pipeline.sh --cores 32 --resources mem_gb=128 -- timetree
```

`timetree` also prepares missing calibrations. Set `phylogeny.dating.enabled: true`
to include dating in `phylogeny`. TimeTree responses are cached across runs;
to refresh, archive `resources/timetree_cache/` and rerun with a new `run_name`.

## Manual calibrations

Create `input/calibrations.tsv` and set `phylogeny.dating.calibration_source: file`
in [config/config.yaml](../config/config.yaml), then run `timetree`.
Each row identifies an ancestor by at least two exact tip labels, separated by
commas. Their most recent common ancestor receives the specified age bounds.
Replace this **example** with suitable species, ages, and a source citation:

```text
taxa	min_age_ma	max_age_ma	source
Species_A,Species_B	90	110	Citation or calibration record
```

Bounds must be positive, with `min_age_ma <= max_age_ma`; equal bounds fix an age.
Unknown tips, duplicate ancestors, and conflicting bounds are rejected.
Calibrations must apply to every requested species set.

## Check the results

Outputs are under `results/<run_name>/phylogeny/<set>/dating/`:

| File | Use |
| --- | --- |
| `species_tree.dated.nwk`, `node_ages.tsv` | Dated tree and internal node ages in Ma |
| `provenance.json` | Settings, chosen smoothing, `review_status`, and `diagnostics` |
| `calibrations.resolved.tsv` | Applied bounds and corresponding node labels |
| `cross_validation.tsv` | Scores from treePL's native cross-validation |
| `treepl.prime.config.txt`, `treepl.config.txt` | Executed options for optimizer selection and dating |

`checks_passed` means numerical checks passed. For `needs_review`, read
`diagnostics`: extend the smoothing grid if its optimum is at a boundary.
If the job fails, read
`logs/<run_name>/phylogeny/<set>/dating.log`; native logs remain in `dating/treepl_runs/`.

Ages are point estimates conditional on the calibrations and model; gene-tree
and calibration uncertainty is not propagated. Native time-branch values are preserved.

## Optional overrides

Smoothing controls the penalty on differences in branch-specific rates. The
workflow runs `prime` once to obtain optimizer settings, then runs treePL with
`cv` once. treePL selects smoothing by native leave-one-out cross-validation and
performs the final fit in the same execution.

Normally omit `phylogeny.dating.treepl`. Add only settings you need to change:

| Setting | When to change it |
| --- | --- |
| `smooth` | Set a reviewed positive value to skip CV; `null` selects by CV. |
| `cvstart`, `cvstop`, `cvmultstep` | Adjust the native grid: `1000.0` to `0.1`, multiplying by `0.1`. |
| `lfiter`, `pliter`, `cviter` | Native optimization iteration counts, defaulting to `3`, `5`, and `3` in the pinned source. |
| `thorough` | Enable extended native optimization; off by default and potentially slow. |

These names map directly to treePL options. The supported settings are defined
in the [dating wrapper](../workflow/scripts/date_phylogeny.py).
Dating-only changes reuse species-tree inference; use a new `run_name` to retain
separate analyses.

## Resources

`date_busco_species_tree` requires **one CPU** and requests **4 GB per species set**.
CV can be expensive for large trees. Check runtime and memory in
`logs/<run_name>/phylogeny/<set>/benchmarks/dating.tsv` and adjust memory using
[rule overrides](running.md#resource-budgets). Keep the thread count at one.

## Method references

For Methods, report the treePL revision, calibration sources, `prime` optimizer
settings, CV grid or fixed `smooth`, site count, and seed. The saved configuration
files and `provenance.json` record these values; `numsites` is the sum of retained
trimmed alignment lengths.

[treePL options](https://github.com/blackrim/treePL/wiki/Run-Options) ·
[treePL method](https://doi.org/10.1093/bioinformatics/bts492) ·
[Coalescent branch lengths and dating](https://doi.org/10.1093/sysbio/syag038)
