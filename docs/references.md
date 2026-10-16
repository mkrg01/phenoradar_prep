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

Set `odb.ncbi_tax_id` to a mapping level containing all build species in the
selected release. Narrower clades define finer OGs. BUSCO lineage is independent.

| NCBI Taxonomy ID | Clade |
| --- | --- |
| `33090` | Viridiplantae |
| `3193` (default) | Embryophyta |
| `4447` | Liliopsida |
| `38820` | Poales |
| `71240` | Eudicots |

Check the [OrthoDB v12 tree](https://data.orthodb.org/v12/tree), or the tree for
your selected release. A clade change requires a new build and new mappings.

### Preparing and verifying the reference

Builds prepare references under `resources/orthodb/{version}_{ncbi_tax_id}/`.
An existing reference is kept unchanged and reused for that release and clade.

For advance setup and checksum verification, use a saved build's config and
its reference path:

```bash
sbatch --cpus-per-task=1 --mem=40G run_pipeline.sh --configfile results/leaf/pipeline.yaml -- references
python workflow/scripts/verify_odb_reference.py \
  --reference resources/orthodb/v12_3193/reference.json
```

Mapping requires network access and an ODB-mapper supporting the selected
release; fallback to another release is rejected.

New mapping jobs and retries share a pinned reference checksum. Reference
checksums are recorded internally; missing historical records remain unknown.

To change releases, edit `odb.version` and
[create a new build](datasets.md#update-orthodb-for-all-samples). Matching results
are reused; other samples are remapped automatically.

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
