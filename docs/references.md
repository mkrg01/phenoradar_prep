# Reference data

[Documentation](index.md)

Missing references are prepared automatically under `resources/`. Completed
snapshots are reused without automatic updates; keep them with analysis records.

| Stage | Reference |
| --- | --- |
| Metadata preparation | NCBI taxonomy |
| `references`, ODB mapping | OrthoDB v12 at `odb.node` |
| `kegg_references`, KEGG analysis | KOfam/KEGG |
| TimeTree retrieval | [Cached API responses](dating.md#timetree-calibrations) |

## Taxonomy reference

If `resources/taxonomy/taxa.sqlite` is missing, the workflow downloads NCBI taxonomy
and builds an ETE4-compatible SQLite snapshot. To bootstrap from a local database:

```yaml
taxonomy:
  source: /path/to/existing/taxa.sqlite
```

This copies the source only when the snapshot is missing. A configured missing
source fails. The `.json` sidecar records provenance and checksum.

For manual offline preparation:

```bash
python workflow/scripts/snapshot_taxonomy.py \
  --source /path/to/existing/taxa.sqlite \
  --destination resources/taxonomy/taxa.sqlite
```

To refresh taxonomy, archive `resources/taxonomy/` and rerun with a new `run_name`.
Changing `taxonomy.source` alone does not replace an existing snapshot.

## OrthoDB reference

### Choosing an OrthoDB node

`odb.node` is the NCBI Taxonomy ID of an OrthoDB v12 mapping level. Choose a
supported clade containing all species in the dataset; a narrower level defines
finer OGs. The default `3193` covers land plants. It is independent of the BUSCO
`phylogeny.lineage` setting.

| `odb.node` | Clade |
| --- | --- |
| `33090` | Viridiplantae (green plants) |
| `3193` | Embryophyta (land plants) |
| `4447` | Liliopsida (monocots) |
| `38820` | Poales |
| `71240` | Eudicots |

Look up clades in [NCBI Taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy) and the
[OrthoDB v12 tree](https://data.orthodb.org/v12/tree). Not every NCBI ID has a
mapping reference. Confirm supported nodes in an [ODB-mapper environment](../workflow/envs/odb.yaml):

```bash
(
  export ODBAPI_URL_VERSION=v12
  export ODBMAPPER_WORK="$PWD/work/odb-node-lookup"
  ODB-mapper SETUP
  ODB-mapper DOWNLOAD '?'
)
```

`'?'` lists nodes without downloading sequences; `'?plants'` restricts the list.
Keep the quotes. See the [OrthoDB v12 guide](https://www.ezlab.org/orthodb_v12_userguide.html)
for OG definitions. Changing `odb.node` requires new mapping/expression results;
use a new `run_name` to retain the previous analysis.

### Preparing and verifying the reference

Snapshots live in `resources/orthodb/v12_<node>/`. To prepare one separately,
after adjusting your [Slurm settings](running.md#slurm):

```bash
sbatch --cpus-per-task=1 --mem=40G \
  run_pipeline.sh --configfile config/mydata.yaml -- references
```

ODB-mapper requires network access during mapping too. Default free-space checks
are 200 GiB for reference preparation and 750 GiB for mapping; these configurable
thresholds are not measured requirements. Nonlocal work storage requires
`odb.allow_nonlocal: true`.

To verify every reference checksum, replace `3193` with your node:

```bash
python workflow/scripts/verify_odb_reference.py \
  --reference resources/orthodb/v12_3193/reference.json
```

For a deliberate refresh, archive the node's snapshot and use a new `run_name`.

## KOfam and KEGG reference

Initial setup downloads profiles and `ko_list` from
[KOfam](https://www.genome.jp/ftp/db/kofam/) plus KO-to-MODULE/PATHWAY maps from
KEGG REST. To prepare without assemblies or annotation:

```bash
./run_pipeline.sh --cores 1 --resources mem_gb=4 -- kegg_references
```

The snapshot is `resources/kegg/snapshot_v1/`; verified downloads are cached in
`resources/kegg/downloads/`. Retries reuse completed downloads. Existing snapshots
need no network access and must not be modified in place.

For local profiles and a matching `ko_list`, with the destination absent:

```bash
python workflow/scripts/prepare_kegg_reference.py \
  --profiles-dir /path/to/kofam/profiles --ko-list /path/to/kofam/ko_list \
  --reference-dir resources/kegg/snapshot_v1 \
  --release YOUR_KOFAM_RELEASE_OR_DOWNLOAD_DATE
```

For fully offline setup, also supply `--module-links` and `--pathway-links` with
headerless two-column responses from `https://rest.kegg.jp/link/module/ko` and
`https://rest.kegg.jp/link/pathway/ko`. Otherwise the helper retrieves those maps.

```bash
python workflow/scripts/verify_kegg_reference.py \
  --reference resources/kegg/snapshot_v1/reference.json
```

To refresh, archive all of `resources/kegg/`, including the download cache, and
rerun with a new `run_name`. Keeping the cache reuses the old downloaded data.
