# Reference data

[Back to README](../README.md)

`resources/` holds reusable reference snapshots, separate from the input assemblies and analysis outputs.
Provide a taxonomy snapshot before running `prepare`; the workflow prepares OrthoDB automatically
when it is needed. Neither database is distributed with the repository.

Run all commands from the repository root. For the batch command below, activate
the workflow environment and create `logs/` before submission, as described in
the [Slurm instructions](running.md#slurm-run-the-workflow-in-one-allocation).

## Taxonomy reference

Provide an existing ETE4-compatible `taxa.sqlite` database. If you do not have one,
first create it with ETE4 in an environment containing that package; the initial
taxonomy download is a separate setup step.

Create a local snapshot, replacing the source path with your database:

```bash
python workflow/scripts/snapshot_taxonomy.py \
  --source /path/to/existing/taxa.sqlite \
  --destination resources/taxonomy/taxa.sqlite
```

This command uses SQLite's backup API, records a checksum, and refuses to
overwrite an existing snapshot. It does not require ETE4 in the environment
running the snapshot command. Metadata processing uses the snapshot without
downloading or updating taxonomy.

To update taxonomy, create a new snapshot and change `taxonomy.database`.

## OrthoDB reference

The full workflow downloads and prepares OrthoDB automatically if the configured
reference snapshot does not exist. To prepare it separately:

```bash
sbatch --partition=YOUR_PARTITION --cpus-per-task=1 --mem=40G \
  run_pipeline.sh --configfile config/mydata.yaml -- references
```

ODB-mapper contacts the OrthoDB API at startup, so mapping also requires network
access after reference preparation. Reference preparation checks for at least
200 GiB of free disk space; mapping checks for at least 750 GiB. These are
configurable thresholds, not measured requirements for every dataset. ODB work
uses local storage by default; using another filesystem requires
`odb.allow_nonlocal: true`. At most one ODB reference/mapping step runs at a time.

To update OrthoDB, use a new reference directory and change `odb.reference_dir`.
Completed reference snapshots are not automatically updated. Mapping checks the
reference inventory checksum, file presence, and file sizes. To verify all
reference file checksums:

```bash
python workflow/scripts/verify_odb_reference.py \
  --reference resources/orthodb/v12_3193/reference.json
```

Adjust the path when using another reference directory.
