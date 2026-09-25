# Reference data

[Documentation](index.md)

Missing references are prepared automatically under `resources/`. Completed
snapshots are reused without updates across run names; keep them with analysis
records. TimeTree uses a separate [response cache](dating.md#timetree-calibrations).

## Taxonomy reference

Metadata preparation builds `resources/taxonomy/taxa.sqlite` from NCBI taxonomy
when absent. Its JSON sidecar records provenance/checksum. To refresh, archive
`resources/taxonomy/` and rerun with a new `run_name`.

## OrthoDB reference

### Choosing an OrthoDB node

`odb.node` is an NCBI Taxonomy ID supported as an OrthoDB v12 mapping level.
Choose a clade containing all dataset species; narrower levels define finer OGs.
It is independent of BUSCO's `phylogeny.lineage`.

| `odb.node` | Clade |
| --- | --- |
| `33090` | Viridiplantae |
| `3193` (default) | Embryophyta |
| `4447` | Liliopsida |
| `38820` | Poales |
| `71240` | Eudicots |

Not every NCBI ID is supported. Check the
[OrthoDB tree](https://data.orthodb.org/v12/tree), or list nodes inside an
[ODB-mapper environment](../workflow/envs/odb.yaml):

```bash
(
  export ODBAPI_URL_VERSION=v12
  export ODBMAPPER_WORK="$PWD/work/odb-node-lookup"
  ODB-mapper SETUP
  ODB-mapper DOWNLOAD '?'
)
```

Keep the quotes: `'?'` lists nodes without sequence downloads; `'?plants'`
restricts the list. Changing nodes requires new mapping/expression results;
use a new `run_name` to preserve earlier analyses.

### Preparing and verifying the reference

Snapshots live in `resources/orthodb/v12_<node>/`. To prepare separately:

```bash
sbatch --cpus-per-task=1 --mem=40G run_pipeline.sh --configfile builds/build001/pipeline.yaml -- references
```

Mapping also needs network access. Verify checksums with your node:

```bash
python workflow/scripts/verify_odb_reference.py \
  --reference resources/orthodb/v12_3193/reference.json
```

To refresh, archive the node's snapshot and use a new `run_name`.

### Reusing existing ODB results

Set `odb.existing_results` to a snapshot directory, relative to the repository
root (or an absolute path). Keep it inside the repository for container runs:

```yaml
odb:
  existing_results: resources/odb_existing/tlight
  node: 3193
```

The directory must contain `annotations.tsv` and `snapshot.json`. Schema version 1
records `version` (`v12`), `node`, `proteins` (one record per species with `species`,
`odb_species`, and the input FASTA's `sha256`), and `annotations` (with
`path: annotations.tsv` and `sha256`). Existing recorded absolute input paths are
provenance only; the original files need not remain at those paths.

For mixed existing/new species, use [incremental mapping](datasets.md#incremental-odb-outside-the-dataset-interface).
The following describes the default strict import mode (`odb.incremental: false`).

Use the normal `mapping` or `all` target. Both bypass reference downloads and
ODB-mapper when this setting is supplied. `null` selects new mapping. The explicit
`references` target still prepares a reference and is unnecessary for imports.

CDS translation still runs: every selected protein FASTA must have the same SHA256
as its recorded original input. The import also checks the OrthoDB version/node,
species coverage, and annotations checksum. Mismatches stop the import; they never
silently trigger remapping. A subset of the recorded species is supported, with
other species' annotations omitted. Output tables retain their normal formats;
`orthogroups/mapping/merge_qc.json` records `mode: existing`, excluded rows, and
snapshot provenance. Abundance-only changes reuse the completed mapping.

## KOfam and KEGG reference

Setup downloads [KOfam](https://www.genome.jp/ftp/db/kofam/) profiles/`ko_list`
and KEGG REST KO-to-MODULE/PATHWAY maps. Prepare without assemblies:

```bash
./run_pipeline.sh --cores 1 --resources mem_gb=4 --configfile analyses/analysis001/pipeline.yaml -- kegg_references
```

The snapshot is `resources/kegg/snapshot_v1/`; retries reuse downloads in
`resources/kegg/downloads/`. Snapshots need no network and must not be edited.

For local profiles and matching `ko_list`, with the destination absent:

```bash
python workflow/scripts/prepare_kegg_reference.py \
  --profiles-dir /path/to/kofam/profiles --ko-list /path/to/kofam/ko_list \
  --reference-dir resources/kegg/snapshot_v1 --release YOUR_RELEASE_OR_DATE
```

For offline setup, also pass `--module-links` and `--pathway-links` as headerless
two-column responses from `https://rest.kegg.jp/link/module/ko` and
`https://rest.kegg.jp/link/pathway/ko`; otherwise those maps are retrieved.

```bash
python workflow/scripts/verify_kegg_reference.py \
  --reference resources/kegg/snapshot_v1/reference.json
```

To refresh, archive all of `resources/kegg/`, including the download cache,
and rerun with a new `run_name`. Keeping the cache reuses old data.
