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
`logs/<analysis>/taxonomy_reference.log`.

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
`analysis` name to retain the results produced with the previous reference.

## OrthoDB reference

The full workflow downloads and prepares OrthoDB automatically if
`resources/orthodb/v12_<node>/reference.json` does not exist. `<node>` comes from
`odb.node`, so changing nodes selects a separate reference directory automatically.
To prepare it separately:

```bash
sbatch --partition=YOUR_PARTITION --cpus-per-task=1 --mem=40G \
  run_pipeline.sh --configfile config/mydata.yaml -- references
```

ODB-mapper contacts the OrthoDB API at startup, so mapping also requires network
access after reference preparation. Reference preparation checks for at least
200 GiB of free disk space; mapping checks for at least 750 GiB. These are
configurable thresholds, not measured requirements for every dataset. ODB work
uses local storage by default; using another filesystem requires
`odb.allow_nonlocal: true`.

To refresh a node, archive its existing reference directory and use a new
`analysis` name before running again. Mapping checks the
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
./run_pipeline.sh --software-deployment-method conda \
  --cores 1 --resources mem_gb=4 -- kegg_references
```

Completed downloads and their URL, retrieval time, byte count, and SHA-256
records stay in `resources/kegg/downloads/`. Retrying after a failed download
reuses verified completed files and restarts the interrupted file. Extraction
and validation must all succeed before the snapshot is published. Existing
snapshots are reused without network access or automatic updates; corrupt
snapshots or cached downloads cause an error rather than being silently replaced.
Preparation logs are in `logs/<analysis>/kegg/reference_prepare.log`.

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
cache, and run again with a new `analysis` name to preserve previous results.
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
