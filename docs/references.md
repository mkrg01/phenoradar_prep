# Reference data

[Documentation](index.md)

Missing references are prepared automatically under `resources/` when needed.
Completed snapshots are reused without automatic updates. Databases are not
distributed with the repository.

| Target or stage | Reference used |
| --- | --- |
| Metadata preparation | NCBI taxonomy snapshot |
| `references`, ODB mapping | OrthoDB v12 at `odb.node` |
| `kegg_references`, KEGG analysis | KOfam/KEGG snapshot |
| TimeTree calibration retrieval | Cached API responses; see [dating](dating.md#timetree-calibrations) |

Run commands from the repository root; see [Slurm setup](running.md#slurm) for submission.

## Taxonomy reference

When `resources/taxonomy/taxa.sqlite` is missing, the workflow downloads NCBI's
`taxdump.tar.gz` and builds an ETE4-compatible SQLite snapshot before selecting
metadata. The first build needs network access to
`https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/`; its log is
`logs/<run_name>/taxonomy_reference.log`.

An existing snapshot is unchanged even if the bootstrap source or workflow code
changes. The `references` target prepares only OrthoDB.

To use an existing local ETE4-compatible SQLite database as the source instead
of downloading taxonomy, add this to your configuration:

```yaml
taxonomy:
  source: /path/to/existing/taxa.sqlite
```

The source is copied with SQLite's backup API only when the destination is
missing. `taxonomy.source: null` (the default) selects the NCBI download. An
explicitly configured missing source is an error; it does not silently switch
to a different taxonomy release.

The database is validated before publication; failed builds leave no partial
snapshot. Rerun the same command to retry. `<database>.json` records
the source, creation time, and database checksum; downloaded builds additionally
record the taxdump checksum and ETE4 version.

Manual offline snapshot creation remains available:

```bash
python workflow/scripts/snapshot_taxonomy.py \
  --source /path/to/existing/taxa.sqlite \
  --destination resources/taxonomy/taxa.sqlite
```

The helper requires no ETE4, records a checksum, and refuses to overwrite an
existing snapshot.

For an intentional taxonomy refresh, archive the existing `resources/taxonomy/`
directory with its analysis records. The workflow creates a new snapshot at the
same fixed path on the next run, using `taxonomy.source` if configured. Use a new
`run_name` to retain the results produced with the previous reference.

## OrthoDB reference

### Choosing an OrthoDB node

`odb.node` specifies the taxonomic level at which OrthoDB orthogroups (OGs) are
defined. Its numeric identifier comes from **NCBI Taxonomy**: the default
[`3193` is Embryophyta (land plants)](https://www.ncbi.nlm.nih.gov/taxonomy/3193).
OrthoDB uses that ID for the mapping level, including the `at3193` suffix in OG
identifiers. See the [OrthoDB v12 guide](https://www.ezlab.org/orthodb_v12_userguide.html)
for the definitions of levels of orthology and OG identifiers.

The workflow passes this value to both `ODB-mapper DOWNLOAD <node>` when
preparing the reference and `ODB-mapper MAP <label> <manifest> <node>` when
mapping proteins. Set it explicitly for your dataset. The separate
`phylogeny.lineage` setting (default `embryophyta_odb12`) describes the BUSCO
tables used for phylogeny; changing it does not change `odb.node`.

To choose another value:

1. Look up the clade name and numeric Taxonomy ID in
   [NCBI Taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy).
2. Confirm that the ID is a supported **mapping level in OrthoDB v12**.
   The [v12 API tree](https://data.orthodb.org/v12/tree) gives clade names (`name`),
   IDs (`key`), and their hierarchy (`children`). Cross-check against the
   downloadable nodes reported by your installed ODB-mapper, as shown below;
   an arbitrary NCBI Taxonomy ID may have no mapping reference.
3. Choose a supported clade that contains all species being compared. A broader
   level defines more inclusive OGs; a narrower level gives finer groupings.
   For a dataset spanning monocots and eudicots, `3193` covers both.

Run this in an environment containing ODB-mapper, as defined in
[workflow/envs/odb.yaml](../workflow/envs/odb.yaml):

```bash
(
  export ODBAPI_URL_VERSION=v12
  export ODBMAPPER_WORK="$PWD/work/odb-node-lookup"
  ODB-mapper SETUP
  ODB-mapper DOWNLOAD '?'
)
```

Here `DOWNLOAD '?'` lists supported node IDs without downloading reference
sequences. Replace `'?'` with `'?plants'` to list only plant nodes. Keep the
quotes so the shell passes the question mark literally. The lookup initializes
a small work directory and requires network access to the OrthoDB API.
`ODB-mapper HELP++` documents the query syntax. In the pinned OrthoLoger 3.8.1
package, the supported IDs come from `etc/v12/orthomapper_levels.sh`.

Examples verified against the v12 API tree and OrthoLoger 3.8.1 node definitions:

| `odb.node` | Clade | Scope |
| --- | --- | --- |
| `33090` | Viridiplantae | Green plants, including green algae and land plants |
| `3193` | Embryophyta | Land plants; the default |
| `4447` | Liliopsida | Monocots |
| `38820` | Poales | The order containing grasses |
| `71240` | eudicotyledons | Eudicots |

For example, for an analysis restricted to eudicots:

```yaml
run_name: eudicots
odb:
  node: 71240
```

When changing the node for an existing dataset, also use a new `run_name`
to retain earlier results. The node changes the OG definitions and mapping
reference, so regenerate the mapping and downstream expression tables.

### Preparing and verifying the reference

The full workflow downloads and prepares OrthoDB automatically if
`resources/orthodb/v12_<node>/reference.json` does not exist. `<node>` comes from
`odb.node`, so changing nodes selects a separate reference directory automatically.
To prepare it separately with Slurm, adjust the
[allocation settings](running.md#slurm) for your computing environment:

```bash
sbatch --cpus-per-task=1 --mem=40G \
  run_pipeline.sh --configfile config/mydata.yaml -- references
```

ODB-mapper contacts the OrthoDB API at startup, so mapping also requires network
access after reference preparation. Reference preparation checks for at least
200 GiB of free disk space; mapping checks for at least 750 GiB. These are
configurable thresholds, not measured requirements for every dataset. ODB work
uses local storage by default; using another filesystem requires
`odb.allow_nonlocal: true`.

To refresh a node, archive its existing reference directory and use a new
`run_name` before running again. Mapping checks the
reference inventory checksum, file presence, and file sizes. To verify all
reference file checksums:

```bash
python workflow/scripts/verify_odb_reference.py \
  --reference resources/orthodb/v12_3193/reference.json
```

Replace `3193` with your configured node when verifying another reference.

## KOfam and KEGG reference

The workflow automatically creates `resources/kegg/snapshot_v1/` when a KEGG
target first needs it. It downloads `profiles.tar.gz` and `ko_list.gz` from the
[official KOfam distribution](https://www.genome.jp/ftp/db/kofam/) and KO-to-MODULE
and KO-to-PATHWAY tables from KEGG REST. Initial setup needs network access and
space for the archives, extracted profiles, and final snapshot.

To prepare the reference separately, without assemblies or annotation:

```bash
./run_pipeline.sh \
  --cores 1 --resources mem_gb=4 -- kegg_references
```

Completed downloads and their URL, retrieval time, byte count, and SHA-256
records stay in `resources/kegg/downloads/`. Retrying after a failed download
reuses verified completed files and restarts the interrupted file. Extraction
and validation must all succeed before the snapshot is published. Existing
snapshots are reused without network access or automatic updates; corrupt
snapshots or cached downloads cause an error rather than being silently replaced.
Preparation logs are in `logs/<run_name>/kegg/reference_prepare.log`.

For manual setup using locally extracted profiles and a matching `ko_list`, the
offline helper remains available. The fixed destination must be absent:

```bash
python workflow/scripts/prepare_kegg_reference.py \
  --profiles-dir /path/to/extracted/kofam/profiles \
  --ko-list /path/to/extracted/kofam/ko_list \
  --reference-dir resources/kegg/snapshot_v1 \
  --release YOUR_KOFAM_RELEASE_OR_DOWNLOAD_DATE
```

This manual command fetches the small KO-to-MODULE and KO-to-PATHWAY link tables from KEGG REST.
For offline setup, also pass `--module-links /path/to/ko_module_links.tsv` and
`--pathway-links /path/to/ko_pathway_links.tsv`. These are headerless two-column
responses from `https://rest.kegg.jp/link/module/ko` and
`https://rest.kegg.jp/link/pathway/ko`, respectively.

Either link direction is accepted. Duplicate memberships are removed;
`path:koNNNNN` and `path:mapNNNNN` normalize to `mapNNNNN`. Multiple memberships
are retained. The complete supplied membership tables are saved, including KOs
without a searched profile.

The snapshot contains `profiles/`, `ko_list`, raw mappings, normalized
`ko_modules.tsv` and `ko_pathways.tsv`, `files.json`, and `reference.json`.
Checksums, source/retrieval information, profile counts, and the release label
are recorded. Existing snapshot directories are never overwritten. To refresh
deliberately, archive the whole `resources/kegg/` directory, including its download
cache, and run again with a new `run_name` to preserve previous results.
Keeping the old download cache would reuse the old downloaded data.

Snapshots must remain immutable. The workflow verifies all file checksums once
before annotation. Species jobs also check the inventory hash, file set, and
sizes, without repeatedly hashing all HMMs. Explicit full verification is:

```bash
python workflow/scripts/verify_kegg_reference.py \
  --reference resources/kegg/snapshot_v1/reference.json
```

Changing files inside an existing snapshot is unsupported, including changes
that preserve timestamps. Keep each completed snapshot with its analysis records.
