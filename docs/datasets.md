# Reusable builds and independent analyses

[Documentation](index.md) · [Configuration](configuration.md)

`run_build.sh` produces reusable species artifacts through ODB mapping.
`run_analysis.sh` selects species from a completed build and produces expression,
alignment, phylogeny, and PhenoRadar outputs. Build and analysis names are independent; analysis names starting with `build_`
are reserved for build output directories.

## Prepare a manually curated dataset

Maintain `input/metadata.tsv` yourself. Include **all desired species**, one run
per species, whether the reads come from NCBI or local files. A manually curated
[run exclusion list](#manually-excluding-unusable-accessions) can omit known
unusable runs without deleting their metadata rows. The workflow does not search
for runs or infer biological inclusion decisions.

Required columns are `scientific_name`, `run`, and positive NCBI `taxid`.
Use AMALGKIT-compatible run metadata; additional columns pass through to
GeneGalleon. Species/run IDs must be unique after normalization. Species IDs
replace spaces with underscores; avoid collisions between hyphens and underscores.

For local reads also supply `private_file: yes`, `lib_layout: single` or `paired`,
`read1_path`, and, for paired reads, `read2_path`. These are TSV column values,
not YAML. Relative FASTQ paths are resolved against the metadata file's directory.
Paths are frozen and reads are made available inside the GeneGalleon workspace.
Use a new run ID when the read content changes.

Edit `config/build.yaml` and `config/analysis.yaml` directly. The build and
analysis commands load these files by default; `prepare` freezes the settings
for each run. Paths are relative to the repository root. `WORKFLOW_PYTHON`
(or the older `DATASET_PYTHON`) can select the host Python with PyYAML.
Install the [workflow environment](../environment.yaml) and prepare the
[workflow container](containers.md) before submitting mapping or analysis jobs.

## Pinned GeneGalleon source and SIF

`build.yaml` pins `genegalleon.version`, the exact source `revision`, and an
OCI image digest in `image_uri`. Missing source/SIF files are fetched when a
build needs assembly, BUSCO, or quantification. A reuse-only build does not fetch
them. `plan` never downloads software. To fetch ahead of time:

```bash
./run_build.sh fetch-software
```

Downloads are cached under `genegalleon.cache_dir`. Source and image checksums
are frozen with the build. An HTTPS SIF URL requires `image_sha256`; OCI pulls
record the produced SIF hash. `repository` and `image` allow explicit local
overrides; an existing `<repository>/genegalleon.sif` is also supported.
Update the version, revision, and image digest together for a deliberate upgrade.
Native products record their conditions; incompatible reuse is reported as a
conflict. Use a separate `store` for a deliberate rebuild with different conditions.
Legacy imports retain their known provenance and are not retroactively assigned
the current software version.

## Build through mapping

```bash
./run_build.sh plan
./run_build.sh prepare --name expansion001
./run_build.sh submit --build builds/expansion001 --until busco --dry-run
./run_build.sh submit --build builds/expansion001 --until busco
./run_build.sh status --build builds/expansion001

# After inspection, schedule any remaining prerequisites and mapping:
./run_build.sh submit --build builds/expansion001 --until mapping
```

Endpoints are `assembly`, `busco`, `quant`, and `mapping` (default). Each runs
missing prerequisites. Assembly includes longest-CDS generation. BUSCO always
retains its full table, regardless of future tree settings. No BUSCO threshold,
trait, or species-list filtering applies to a build: every species remaining after explicit run exclusions must
finish all required products before the build is complete.

The plan reports `reuse`, `pending`, or `conflict` for species products. Mapping
reuse requires translated protein hashes: before translation its state is
`check_after_translation` (or `pending_inputs` without CDS). The subsequent
incremental ODB plan records `planned_reuse`/`planned_mapping`; completion gives
`reuse`. An unresolved conflict stops preparation/submission. Reports include
BUSCO completeness once available, for manual inspection.

`builds/<id>/` holds frozen metadata/configuration, the persistent GeneGalleon
workspace, and job records. Mapping results live in `results/build_<id>/`.
After successful mapping, validation checks full metadata coverage, CDS/protein
identity, BUSCO/quant files, and mapping gene coverage. Only then is
`builds/<id>/completed.json` published. Analysis rejects a missing or changed
completion record or product. If an independently run mapping finished before
its completion record was written, `run_build.sh complete --build builds/<id>`
validates and publishes it without running analyses.

For a pilot, put species IDs in a file and add `--species-list pilot.txt` to
`submit`. Only those species' missing upstream tasks are submitted; no mapping
controller or completion record is produced for an incomplete pilot. Later
submit without the pilot list to finish the build.

## Run an analysis

Set `build` in `analysis.yaml`, or override it with `--build` when preparing.
Set BUSCO acceptance, species selection, traits, and optional analyses here:

```bash
./run_analysis.sh plan --build builds/expansion001
./run_analysis.sh prepare --build builds/expansion001 --name carnivory001
./run_analysis.sh submit --analysis analyses/carnivory001 --dry-run
./run_analysis.sh submit --analysis analyses/carnivory001
./run_analysis.sh status --analysis analyses/carnivory001
```

Analysis uses verified build proteins and a selected-species subset of the
completed mapping database. Its DAG contains no assembly, BUSCO, `translate_cds`,
or ODB-mapper producer. Missing build products cause an error; they are not rebuilt.
KO annotation remains an optional analysis branch.

`inputs.species_trait` names the trait TSV (use `null` when unused).
`selection.species_list: true` enables the file named by `inputs.species_list`.
`exclude_species` removes exact species IDs before downstream computation.
`inputs.calibrations` supplies optional manual dating constraints.
Use a new analysis name after changing any of these inputs or scientific settings.
Multiple analyses can share a completed build; they cannot change its lineage,
genetic code, ODB node, or membership. They inherit those values.

`analyses/<id>/` holds frozen inputs and resolved settings. Results, work, and logs
use the analysis ID under `results/`, `work/`, and `logs/`. The default target
`all` also collects `results/<id>/phenoradar_inputs/` after successful execution.
Use `submit --target phylogeny` or another [analysis target](running.md#targets)
for partial downstream execution. `run --local --cores 8 --mem-mb 32000` runs a
prepared analysis directly with the same container launcher.

## Slurm and retries

Assembly, BUSCO, and quantification run as species arrays. `slurm.concurrency`
limits concurrent species tasks; `slurm.array_size` bounds the largest local
array index (default 1000; must be below the site's `MaxArraySize`). Large arrays
are split into dependent batches with stable species indices. Phases and batches
are connected by `afterok`; mapping starts only when upstream jobs succeed.

Mapping and analysis each have a small controller allocation. Snakemake submits
individual rules as separate Slurm jobs with their own time/memory limits.
`slurm.jobs` controls their concurrency; `slurm.stages.controller` controls the
controller, and `slurm.rules` controls rule resources. Update these entries
in `config/build.yaml`, for example:

```yaml
slurm:
  stages:
    assembly: {cpus: 8, mem_mb: 256000, time: "7-00:00:00"}
  rules:
    odb_map: {cpus: 16, mem_mb: 192000, runtime: 4320}
```

After editing the resource settings, apply them to a retry explicitly:

```bash
./run_build.sh submit --build builds/expansion001 --until mapping \
  --resources config/build.yaml
./run_analysis.sh submit --analysis analyses/carnivory001 \
  --resources config/analysis.yaml
```

Only the current `slurm` section is read from the specified file. Scientific settings and
metadata remain frozen. Each submission stores its resolved resources and
scripts separately, including dry-runs. Retrying submits unfinished species
steps and reuses completed mapping chunks. Active or ambiguous recorded Slurm
submissions block overlapping submissions; inspect/cancel them before retrying.
A controller timeout can leave its independently submitted workers active, so
also inspect those worker jobs before resubmitting. Job logs are under the
build/analysis `jobs/logs/` directory.

## Failure and recovery behavior

Failure never automatically excludes an accession or publishes a completed
build. A successful, validated species/stage is registered for reuse. A failed
species/stage stays `pending` in `status`, with its last attempt recorded as
`failed` under `jobs/status/`; the record includes run accession, species, stage,
job ID, and error. Check `jobs/logs/` and the retained GeneGalleon logs for the
underlying download/assembler error. Slurm timeouts or hard kills may leave the
last-attempt record at `running`; that record alone is not proof of completion.

Submit the **same build** again after inspecting the queue:

```bash
./run_build.sh status --build builds/expansion001
./run_build.sh submit --build builds/expansion001 --until mapping
```

Completed species/stages are skipped. Failed/missing ones are retried, and
unfinished published outputs are moved to `jobs/incomplete/` before retrying.
FASTQ/download and temporary work directories are retained. Missing downstream
prerequisites are included in the new submission. Other species already running
in the same array can finish, but the current `afterok` dependencies stop later
batches/phases when any predecessor fails. Remaining queued jobs must finish or
be cancelled before resubmission; retries are explicit, not an unlimited loop.

**Stage retry and within-assembler checkpoint resume are different.** The pinned
GeneGalleon 0.7.77 uses AMALGKIT `getfastq --redo no` and has recovery code for
retained download/run state; reusable reads are kept when valid. It does not
promise recovery of every partially downloaded byte. Its rnaSPAdes path removes
`rnaspades_output` before starting rnaSPAdes again, so an interrupted assembly is
recomputed for that species. Other assembler behavior is tool-dependent.
For repeated timeouts or memory errors, first consider a resource override;
resubmitting unchanged limits may fail again.

## Manually excluding unusable accessions

Edit [config/excluded_accessions.tsv](../config/excluded_accessions.tsv) manually.
It starts with a header and no exclusions. Entries match the metadata **`run`**
column exactly (SRR/ERR/DRR accessions or your local run IDs):

```tsv
accession	reason
SRR123456	download_failed_repeatedly
ERR987654	assembly_failed_after_review
LOCAL_BAD1	unusable_reads
```

These are illustrative IDs, not built-in exclusions. `accession` is required;
`reason` is optional but recommended. Additional columns such as `stage`,
`notes`, or `reviewed_on` are retained. Duplicate/malformed entries are errors.
Matching is case-sensitive and uses exact IDs, not prefixes or regular expressions.
BioProject/experiment/species identifiers do not exclude their associated runs.

The build config selects the list:

```yaml
excluded_accessions: config/excluded_accessions.tsv
```

Use `null` to disable it. Older configs without this key retain their
old behavior; add the key to enable exclusions. `plan` and `status` show excluded
runs with their reasons. `prepare` removes them before resolving cached products
or requiring FASTQs, and they receive no worker-array indices. `register` also
skips them. Excluding the only run for a species removes that species from the
new build and all analyses based on it; no replacement run is chosen automatically.
The remaining metadata still requires one run per species.

After editing the list, prepare a **new build ID** using the same metadata and
store. The new build reuses eligible completed products. An existing build keeps
its frozen membership; source-list edits neither modify it nor cancel its jobs.
Remove an entry to permit the run in a later build again. Products and historical
outputs are retained. The list selects metadata runs; it does not retroactively
invalidate a CDS reference originally assembled from another run of the same
species. Reference rejection/replacement remains an explicit choice of
`reference_id` or a separate store.

Each build preserves the original `source_metadata.tsv`, the exact exclusion TSV,
its effective `metadata.tsv`, and an `excluded_runs.tsv` audit with species/run/reason.
These records are checksummed with the build. A list may contain accessions absent
from current metadata: keep those decisions so future metadata preparation can
avoid them too. The independent `accession_exclusions.read_exclusions()` reader
can be reused by a future metadata generator. No automatic failure classification
or automatic editing of the exclusion list is performed.

## Updating species and importing existing work

Edit metadata and prepare a **new build ID**. Added species run missing work;
removed species are absent from the new build and its analyses. Existing build
membership and historical results are unchanged. Species receipts in
`resources/dataset_assets/` and mappings in `resources/odb_cache/` persist, so
re-adding a species can reuse them. A changed run can reuse CDS/BUSCO/mapping
and request quantification for the new run.

To register already assembled/quantified species using the [input layout](inputs.md):

```bash
./run_build.sh register \
  --input-dir input --metadata input/metadata.tsv
```

Registration validates artifacts without recomputation. Missing CDS is reported;
missing full BUSCO tables or quantification remain pending. A summary alone
cannot complete the new build. Preserve original data/workspaces referenced by
receipts. If multiple CDS references exist, set the desired `reference_id` in
metadata explicitly.

For old ODB results, set `odb.existing_results` in **build.yaml** to an imported
snapshot; see [ODB imports](references.md#reusing-existing-odb-results). The build
combines matching old snapshots, native cached chunks, and missing species.
Protein/reference mismatches are conflicts rather than silent reuse.

## Migrating old configurations

`config/config.yaml` and `config/dataset.yaml` are replaced by the two phase
configs. Existing stores, snapshots, `datasets/`, and results are not deleted.
Convert saved old settings with:

```bash
python workflow/scripts/migrate_phase_config.py \
  --legacy-dataset config/dataset.local.yaml --legacy-config saved-config.yaml \
  --build-output config/build.migrated.yaml --analysis-output config/analysis.migrated.yaml
```

Alternatively, pass `--legacy-dataset datasets/<old-id>/dataset.json` to recover
its frozen settings/metadata. The converter retains cache paths and writes only
new configuration files; it refuses to overwrite existing files. Review the
generated settings and apply them to `config/build.yaml` and `config/analysis.yaml`,
then select the new completed build in analysis settings.
Register old input snapshots if needed, then prepare a new build. Old schema-1
dataset jobs should be resumed with their original checkout, not modified in
place. `run_dataset.sh` remains a command-name alias for the new build interface;
its old `--until all` combined execution has been removed.
