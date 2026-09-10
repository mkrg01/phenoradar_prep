# Reference data

[Back to README](../README.md)

`resources/` holds reusable reference snapshots, separate from the input assemblies and analysis outputs.
The workflow prepares missing taxonomy, OrthoDB, and KOfam/KEGG references
automatically when they are needed. These databases are not distributed with
the repository, and their storage locations are fixed by the workflow.

The optional KEGG branch uses its own immutable KOfam/KEGG snapshot at
`resources/kegg/snapshot_v1/`. It is created on first use, independently of ODB.
The `kegg_references` target prepares it separately; `references` prepares only
OrthoDB. See [KEGG setup](kegg.md) for download caching and offline preparation.

Run all commands from the repository root. For the batch command below, activate
the workflow environment and create `logs/` before submission, as described in
the [Slurm instructions](running.md#slurm-run-the-workflow-in-one-allocation).

## Taxonomy reference

When `resources/taxonomy/taxa.sqlite` is missing, the workflow downloads NCBI's
`taxdump.tar.gz` and builds an ETE4-compatible SQLite snapshot before selecting
metadata. This happens automatically for `prepare`, the full workflow, and the
standalone `kegg` target. No extra setup command is required. The first build
needs network access to `https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/`; its log is
`logs/<analysis>/taxonomy_reference.log`.

An existing taxonomy snapshot is used unchanged, even if the bootstrap source
or workflow code changes. It is not downloaded or updated on each run. A dry-run
only schedules creation and does not download anything. The `references` target
continues to prepare only OrthoDB.

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

Downloads and ETE working files are isolated in a temporary directory beside
the destination. The SQLite database is validated before publication, so a failed
download or build leaves no partial database at the fixed path. Resubmit
the same workflow command to retry. Creation also writes `<database>.json` with
the source, creation time, and database checksum; downloaded builds additionally
record the taxdump checksum and ETE4 version.

Manual offline snapshot creation remains available:

```bash
python workflow/scripts/snapshot_taxonomy.py \
  --source /path/to/existing/taxa.sqlite \
  --destination resources/taxonomy/taxa.sqlite
```

This command uses SQLite's backup API, records a checksum, and refuses to
overwrite an existing snapshot. It does not require ETE4 in the environment
running the snapshot command. Metadata processing itself always uses the
prepared snapshot without downloading or updating taxonomy.

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
`odb.allow_nonlocal: true`. CPU and memory budgets determine how many ODB jobs
can run concurrently.

Completed reference snapshots are not automatically updated. To refresh the same
node deliberately, archive its existing reference directory and use a new
`analysis` name before running again. Mapping checks the
reference inventory checksum, file presence, and file sizes. To verify all
reference file checksums:

```bash
python workflow/scripts/verify_odb_reference.py \
  --reference resources/orthodb/v12_3193/reference.json
```

Replace `3193` with your configured node when verifying another reference.
