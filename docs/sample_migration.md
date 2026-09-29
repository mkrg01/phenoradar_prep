# One-time sample ID migration

New builds use `Species_name_accession` for every independently assembled sample.
Biological names and taxids are preserved in separate columns. The runtime adds
these IDs to curated metadata automatically; it does not infer species by splitting
an ID. See [identifiers](inputs.md#identifiers).

`workflow/scripts/migrate_sample_ids.py` converts a legacy registered dataset into
a separate destination. It changes FASTA identifiers, BUSCO sequence identifiers,
abundance `target_id`, and ODB `gene_id` together, preserving sequences, numeric
text, assignments and coordinates. It validates the results and recreates the
CDS/BUSCO/quant registry, translations and ODB cache under their new checksums.
No assembly, BUSCO run, quantification, translation or ODB search is performed.

```bash
python workflow/scripts/migrate_sample_ids.py \
  --metadata input/metadata.tsv \
  --store resources/dataset_assets \
  --products builds/angiosperm_leaf_20260925/products \
  --legacy-results results/run001 \
  --old-cache resources/odb_cache \
  --output migrations/sample_ids_20260928 \
  --workers 8
```

`--limit 2` performs a pilot without publishing a completed migration. Running
without the limit resumes matching validated per-sample receipts. Changed source
identities or conflicting destination contents stop conversion. The source must
have an unambiguous species/run/reference association. Previously excluded but
registered samples are retained; exclusions continue to apply to new builds.

The destination contains `sample_id_map.tsv`, per-sample `receipts/`, transformed
`input/`, `dataset_assets/`, `mapping/`, `odb_cache/`, and `completed.json` only
after the entire migration succeeds. Source artifacts are never edited in place:
they can be hard-linked into historical builds. The destination requires space
for a new copy of the transformed products.

After validation, publish a new build without running scientific tools:

```bash
python workflow/scripts/publish_sample_migration.py \
  --migration migrations/sample_ids_20260928 \
  --name angiosperm_leaf_20260928_samples
```

This prepares a new build using the migrated registry and ODB cache,
materializes the selected inputs, attaches the converted proteins and mapping
tables, and publishes a new completed products bundle. It writes
`build_validation.json` and a proposed `build.yaml` in the migration directory.
Update the default build and analysis settings only after validating this bundle
and confirming reuse of every selected sample. Old build/analysis records remain historical snapshots; they are not edited
or relabelled silently.

The one-time `migrate_sample_ids.py` and `publish_sample_migration.py` tools and
legacy import documentation can be removed after cutover. `sample_identity.py`
and `relabel_sample.py` remain normal build helpers for newly processed samples.
Keep the ID map and validation receipts with the migrated dataset. The destination
also contains the active registered products and caches: retain it while build
settings or registry records reference it, even after removing the migration code. The public
`register` import command may also be removed if no longer needed, while normal
build publication still needs internal registry and cache management.
