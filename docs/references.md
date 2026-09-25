# Reference data

[Documentation](index.md)

Missing references are prepared automatically under `resources/`. Completed
snapshots are shared across builds/analyses without automatic updates. Keep their
provenance with results. TimeTree has a separate [response cache](dating.md#timetree-calibrations).

## Taxonomy reference

Analysis metadata preparation creates `resources/taxonomy/taxa.sqlite` from NCBI
taxonomy when absent, with a provenance/checksum sidecar. To refresh, archive
`resources/taxonomy/` and prepare a new analysis.

## OrthoDB reference

### Choosing an OrthoDB node

`odb.node` in `build.yaml` is an NCBI taxid supported as an OrthoDB v12 mapping level.
Choose a clade containing all build species; narrower levels define finer OGs.
It is independent of `busco.lineage`.

| `odb.node` | Clade |
| --- | --- |
| `33090` | Viridiplantae |
| `3193` (default) | Embryophyta |
| `4447` | Liliopsida |
| `38820` | Poales |
| `71240` | Eudicots |

Check supported nodes in the [OrthoDB tree](https://data.orthodb.org/v12/tree),
or run this in an [ODB-mapper environment](../workflow/envs/odb.yaml):

```bash
(
  export ODBAPI_URL_VERSION=v12
  export ODBMAPPER_WORK="$PWD/work/odb-node-lookup"
  ODB-mapper SETUP
  ODB-mapper DOWNLOAD '?'
)
```

Quoted `'?'` lists nodes without sequences; `'?plants'` restricts the list.
Changing the node requires a new build and new mappings.

### Preparing and verifying the reference

Build prepares `resources/orthodb/v12_<node>/` automatically. For separate setup,
use a prepared build's resolved config:

```bash
sbatch --cpus-per-task=1 --mem=40G run_pipeline.sh --configfile builds/build001/pipeline.yaml -- references
python workflow/scripts/verify_odb_reference.py \
  --reference resources/orthodb/v12_3193/reference.json
```

Mapping also needs network access. To refresh, archive the node's snapshot and
prepare a new build, reviewing the compatibility of retained mapping caches.

### Reusing existing ODB results

Import a snapshot containing `annotations.tsv` and `snapshot.json` once:

```bash
./run_build.sh register --odb-results imports/old_odb --odb-only
```

The snapshot records v12/node, per-species protein SHA256 hashes, and the annotation
checksum. Registration validates it and copies/hard-links it into
`odb.cache_dir/v12_<node>/`, preserving the source. Subsequent builds discover
matching species automatically, map missing ones, and omit excluded species.
Protein mismatches are conflicts; abundance-only changes can reuse mappings.

New builds still translate CDS to verify protein identity. Analysis of a
[completed bundle](datasets.md#copying-a-completed-build-to-another-project)
reuses proteins and mappings directly. `odb.existing_results` remains accepted
for legacy configs but is unnecessary for normal builds.

## KOfam and KEGG reference

KEGG analysis downloads [KOfam](https://www.genome.jp/ftp/db/kofam/) profiles/`ko_list`
and KEGG REST KO-to-MODULE/PATHWAY maps into `resources/kegg/snapshot_v1/`.
Retries reuse `resources/kegg/downloads/`; completed snapshots work offline.
For separate setup, the low-level target is `kegg_references`.

To prepare from matching local profiles/`ko_list` with an absent destination:

```bash
python workflow/scripts/prepare_kegg_reference.py \
  --profiles-dir /path/to/kofam/profiles --ko-list /path/to/kofam/ko_list \
  --reference-dir resources/kegg/snapshot_v1 --release YOUR_RELEASE_OR_DATE
python workflow/scripts/verify_kegg_reference.py \
  --reference resources/kegg/snapshot_v1/reference.json
```

For offline setup, also supply `--module-links` and `--pathway-links`: headerless
two-column responses from `https://rest.kegg.jp/link/module/ko` and
`https://rest.kegg.jp/link/pathway/ko`. Otherwise those maps are downloaded.
To refresh, archive all of `resources/kegg/`, including downloads, and prepare a
new analysis. Keeping the download cache reuses old data.
