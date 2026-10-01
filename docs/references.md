# Reference data

[Documentation](index.md)

Required taxonomy, OrthoDB, and KEGG references are prepared automatically under
`resources/`. Completed snapshots are reused without automatic updates; keep
their provenance with results. GeneGalleon manages upstream references such as
BUSCO. TimeTree uses a separate [response cache](dating.md#timetree-calibrations).

## Taxonomy reference

`resources/taxonomy/taxa.sqlite` is created from NCBI taxonomy when absent.
To refresh, archive `resources/taxonomy/` and prepare a new build or analysis.

## OrthoDB reference

### Choosing an OrthoDB mapping clade

Set `odb.ncbi_tax_id` in `build.yaml` to an OrthoDB v12 mapping level containing
all build species. Narrower clades define finer OGs. This choice is independent
of `busco.lineage`.

| NCBI Taxonomy ID | Clade |
| --- | --- |
| `33090` | Viridiplantae |
| `3193` (default) | Embryophyta |
| `4447` | Liliopsida |
| `38820` | Poales |
| `71240` | Eudicots |

Check the [OrthoDB tree](https://data.orthodb.org/v12/tree) for supported nodes.
Changing the clade requires a new build and new mappings.

### Preparing and verifying the reference

Build prepares `resources/orthodb/v12_<ncbi_tax_id>/` automatically. For advance
setup, use a prepared build's config:

```bash
sbatch --cpus-per-task=1 --mem=40G run_pipeline.sh --configfile results/leaf/pipeline.yaml -- references
python workflow/scripts/verify_odb_reference.py \
  --reference resources/orthodb/v12_3193/reference.json
```

Mapping also needs network access. To refresh, archive the clade's snapshot and
prepare a new build. Reuse only databases with compatible references; see
[database reuse](datasets.md#reusing-completed-databases).

## KOfam and KEGG reference

The KEGG branch downloads KOfam profiles, `ko_list`, and KO-to-MODULE/PATHWAY maps
into `resources/kegg/snapshot_v1/`. Searches use the snapshot's `eukaryote.hal`.
Completed snapshots work offline; download retries reuse `resources/kegg/downloads/`.
The separate setup target is `kegg_references`.

To prepare a new snapshot from local files, include profiles, `ko_list`, and the
original `eukaryote.hal` from the same KOfam release:

```bash
python workflow/scripts/prepare_kegg_reference.py \
  --profiles-dir /path/to/kofam/profiles --ko-list /path/to/kofam/ko_list \
  --reference-dir resources/kegg/snapshot_v1 --release YOUR_RELEASE_OR_DATE
python workflow/scripts/verify_kegg_reference.py \
  --reference resources/kegg/snapshot_v1/reference.json
```

For offline setup, also supply `--module-links` and `--pathway-links` as headerless
two-column responses from KEGG's `/link/module/ko` and `/link/pathway/ko` endpoints.
Otherwise those maps are downloaded.

Snapshots lacking `eukaryote.hal` must be rebuilt from matching inputs. Do not
add a list from a different release to an existing snapshot. To refresh all KEGG
data, archive **all of `resources/kegg/`, including downloads**, and prepare a
new analysis; keeping the download cache reuses old data.
