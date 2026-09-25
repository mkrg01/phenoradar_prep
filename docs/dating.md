# Dating species trees

[Documentation](index.md) · [BUSCO phylogeny](phylogeny.md)

`timetree` uses **treePL** to estimate ages in millions of years (Ma), with
**TimeTree** calibrations by default or manual age bounds. Topology and root stay
fixed; all tips are treated as living species with age zero.

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

## Manual calibrations

Set `inputs.calibrations: input/calibrations.tsv` and
`phylogeny.dating.calibration_source: file` in `analysis.yaml` before preparing.
Each row identifies an ancestor by at least two exact tip labels, comma-separated.
Their most recent common ancestor receives the bounds. Replace this example with
appropriate species, ages, and citations:

```tsv
taxa	min_age_ma	max_age_ma	source
Species_A,Species_B	90	110	Citation or calibration record
```

Bounds must be positive with `min_age_ma <= max_age_ma`; equal bounds fix an age.
Unknown tips, duplicate ancestors, and conflicting bounds are rejected.
Calibrations must apply to every requested species set.

## Check the results

Under `results/<analysis>/phylogeny/<set>/dating/`:

| File | Use |
| --- | --- |
| `species_tree.dated.nwk`, `node_ages.tsv` | Dated tree and ages in Ma |
| `provenance.json` | Settings, smoothing, `review_status`, and `diagnostics` |
| `calibrations.resolved.tsv` | Applied bounds and node labels |
| `cross_validation.tsv` | Native cross-validation scores |
| `treepl.prime.config.txt`, `treepl.config.txt` | Executed optimizer/dating options |

`checks_passed` means numerical checks passed. For `needs_review`, inspect
`diagnostics`; a smoothing optimum at a grid boundary needs review. Failures are
logged in `logs/<analysis>/phylogeny/<set>/dating.log`, with native logs retained
in `dating/treepl_runs/`. Ages are conditional point estimates: gene-tree and
calibration uncertainty is not propagated.

## Optional overrides

The workflow runs treePL `prime` for optimizer settings, then native leave-one-out
cross-validation and final fitting. Normally omit `phylogeny.dating.treepl`;
add only options requiring changes:

| Setting | Use |
| --- | --- |
| `smooth` | Reviewed positive value to skip CV; `null` uses CV |
| `cvstart`, `cvstop`, `cvmultstep` | Default grid: `1000.0` to `0.1`, multiplying by `0.1` |
| `lfiter`, `pliter`, `cviter` | Optimization iteration counts |
| `thorough` | Extended, potentially slow optimization |

See the [wrapper](../workflow/scripts/date_phylogeny.py) for supported options.
Changed scientific settings require a new analysis; tree results are not
shared automatically across analysis IDs.

## Resources

`date_busco_species_tree` requires **one CPU** and defaults to **4 GB per species set**.
Check `logs/<analysis>/phylogeny/<set>/benchmarks/dating.tsv` and adjust memory/time
with [rule overrides](running.md#resource-budgets). Keep one thread.

## Method references

Saved configs/provenance record the treePL revision, sources, optimizer settings,
CV grid or smoothing, retained site count, and seed for reporting Methods.

[treePL options](https://github.com/blackrim/treePL/wiki/Run-Options) ·
[treePL method](https://doi.org/10.1093/bioinformatics/bts492) ·
[Coalescent branch lengths and dating](https://doi.org/10.1093/sysbio/syag038)
