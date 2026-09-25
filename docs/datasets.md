# Incremental RNA-seq datasets

[Documentation](index.md) · [Input contract](inputs.md) · [Slurm execution](running.md)

`run_dataset.sh` integrates manually curated RNA-seq metadata, existing species
products, GeneGalleon, and the downstream workflow. Metadata acquisition and
biological inclusion decisions remain manual. The dataset interface requires
one row and one run per species; it never averages runs or infers sample labels.
The existing prepared-input workflow remains available unchanged.

## Storage and identity

- `resources/software/genegalleon/`: version-pinned source, SIF files, and download receipts, shared across datasets.
- `resources/dataset_assets/`: receipts binding CDS, BUSCO, and quantification to
  a species, CDS SHA256, and run. Legacy files are registered in place; retain
  their original files. New native GeneGalleon outputs remain in their dataset
  workspace, including counts, effective lengths, and tool logs.
- `resources/odb_cache/v12_<node>/`: immutable completed ODB mapping snapshots,
  reusable across datasets, independently of the original mapping batch size.
- `datasets/<name>/`: frozen metadata, configuration, code/image fingerprints,
  stable task indices, GeneGalleon workspace, submission records, and input view.
- `results/<name>/`: results containing only this dataset's selected species.

Never edit registered files in place. A changed registered input is a conflict,
not permission to silently reuse an incompatible output or restart heavy work.
Existing provenance is recorded as legacy; importing does not claim to recover
historical tool versions or prove the original quantification reference beyond
available IDs and retained inputs. Generated products record the GeneGalleon
code/image and effective settings. A metadata-only edit does not invalidate
unchanged species products.

## Configure once

Use the Snakemake environment in [environment.yaml](../environment.yaml), or
install `snakemake-executor-plugin-slurm` alongside your existing Snakemake.
`run_dataset.sh` uses `python`; set `DATASET_PYTHON` to an absolute interpreter
path if that environment is not activated.

Copy [config/dataset.yaml](../config/dataset.yaml) to a local configuration,
then set:

- `genegalleon.version`, `revision`, and `image_uri`: the pinned upstream version,
  full source commit SHA, and matching OCI digest. The defaults select GeneGalleon
  0.7.77; leave `repository` and `image` null for automatic retrieval;
- `analysis_config`: an optional override to the main downstream configuration;
  use `config/reuse_odb.local.yaml` for the local tlight snapshot when available;
- Slurm partition/account, per-stage CPUs/memory/time, concurrency, and the
  downstream job limit. The provided `epyc` and resource values are examples.

The downstream settings in `config/config.yaml` still select BUSCO thresholds,
traits, trees, alignments, and optional analyses. Dataset preparation freezes the
resolved configuration; later edits to the source config apply to later datasets.
Upstream mapping/quantification must use the CDS retained for this reference.
GeneGalleon stage flags, metadata mode, reference selection, and cleanup are
managed by the adapter. Other scalar settings can be set in
`genegalleon.settings`; they are frozen and recorded.

## Pinned GeneGalleon source and SIF

`prepare` automatically fetches missing GeneGalleon dependencies when assembly,
BUSCO, or quantification remains pending. It does not fetch anything for a dataset
that can reuse all upstream products. `plan`, `register`, and `submit --dry-run`
do not download software. Workers use the paths frozen by `prepare`, so every
array task uses the same source and SIF without repeating downloads.

The default pins correspond to the published
[GeneGalleon 0.7.77 source](https://github.com/kfuku52/genegalleon/tree/3a6460e8a8201ab1293db1fec445fd7de063e46e)
and its matching multi-architecture image:

```yaml
genegalleon:
  version: "0.7.77"
  revision: 3a6460e8a8201ab1293db1fec445fd7de063e46e
  image_uri: docker://ghcr.io/kfuku52/genegalleon@sha256:df357fc1857c737df1cdfe17fc12b7e1fcb891353e9bac9e534fb4155b98aff3
  image_sha256: null
  cache_dir: resources/software/genegalleon
  repository: null
  image: null
```

Source archives are fetched by full commit SHA and checked against `VERSION`.
Apptainer (or Singularity) pulls the pinned OCI digest and converts it to a SIF
for the host architecture. This first conversion can require substantial time,
memory, and temporary disk space. Fetch dependencies ahead of dataset preparation
on a suitable node/allocation with:

```bash
./run_dataset.sh fetch-software --config config/dataset.local.yaml
```

This command needs no RNA-seq metadata and submits no analysis jobs. Network
access is needed only for missing dependencies. The source, image, and runtime
blob cache are stored under `genegalleon.cache_dir`, which must be inside the
project. Successfully published caches can be reused offline. Concurrent fetches
are protected by file locks; if another preparation holds the lock, retry after
it finishes. Failed or interrupted downloads never become completed cache entries.

An HTTPS URL for a directly distributed SIF is also supported, with its mandatory
`image_sha256`. OCI tags such as `latest` are rejected: use an immutable digest.
The final SIF checksum, architecture, source commit, source archive checksum,
and available OCI labels are recorded in `dataset.json` under `software_lock`.
Image/source revision and version labels are checked when the image provides them.
OCI-to-SIF conversion may produce different SIF bytes with a different runtime;
each prepared dataset remains bound to its original SIF checksum. Keep its cache
while the dataset is in use. Changed cached files are reported as conflicts.

For an existing checkout/container, set `repository` and `image` to local paths.
If only `repository` is set and `<repository>/genegalleon.sif` exists, it is reused.
These are explicit local overrides; their actual files are fingerprinted rather
than claimed to match the default pins. Without a local SIF, automatic image
retrieval requires the local checkout's `VERSION` to match the configured version.

Update `version`, `revision`, and `image_uri` together when selecting a new release.
Each dependency identity gets a separate cache entry. Existing datasets retain
their frozen paths and fingerprints, and changing software pins does not by itself
rerun completed species products.

## Register completed work

Register existing prepared inputs once, without assembly, BUSCO, or quantification:

```bash
./run_dataset.sh register --config config/dataset.local.yaml \
  --input-dir input --metadata input/metadata.original.tsv
```

Use the metadata describing the original run/reference combinations. Registration
checks CDS IDs, positive finite TPM, equality of abundance/CDS target sets, and
BUSCO counts. Full BUSCO tables additionally check sequence IDs and lineage.
CDS-only and CDS+BUSCO imports are allowed; missing steps remain pending. Species
without CDS are reported as `no_cds`, not marked complete. If the analysis needs
phylogeny, full BUSCO tables must be present or that BUSCO step remains pending.
Files are hashed on registration; later checks use file identity and rehash
changed files. Registration does not require the original FASTQ files.
For newly processed private runs, retained FASTQs are checked against their
recorded hashes when reusing quantification; replacing bytes under the same run
ID is a conflict. Completed outputs can still be reused after raw reads are
removed. Use a new run ID for a different sample.

## Prepare a manually curated dataset

Prepare a TSV containing all desired species, including already processed ones.
It needs `scientific_name`, `run`, and `taxid`, plus the AMALGKIT fields needed for
new work. For local reads specify `private_file=yes`, `lib_layout=single|paired`,
`read1_path`, and (for paired reads) `read2_path`. Relative FASTQ paths resolve
against the metadata file's directory. Biological labels and `is_sampled` /
`exclusion` must already be curated for AMALGKIT. Inclusion in the final dataset
is governed by row membership and configured species/BUSCO filters, not by an
implicit interpretation of AMALGKIT exclusion columns.

```bash
./run_dataset.sh plan --config config/dataset.local.yaml \
  --metadata input/metadata.tsv
./run_dataset.sh prepare --config config/dataset.local.yaml \
  --metadata input/metadata.tsv --name expansion001
```

`plan` reports reuse, pending work, explicit exclusions, and conflicts without
submitting jobs or acquiring metadata. `prepare` freezes the dataset, hashes
new private read inputs, and resolves/fingerprints the pinned GeneGalleon code/image
when upstream work is needed. It does not launch analyses. Use a new dataset name when changing
the species/run selection; existing names and existing result directories cannot
be overwritten. An optional `reference_id` column chooses a registered CDS SHA256
when a species has multiple imported references. No automatic reassembly follows
from adding/changing RNA-seq runs: existing CDS are reused and only missing run
quantification is scheduled. An intentional reference replacement requires a
separately prepared/imported reference and compatible ODB snapshots/cache; an
incompatible old ODB snapshot is reported instead of silently remapping.

## Submit through a chosen endpoint

Inspect generated batch commands before submitting:

```bash
./run_dataset.sh submit --dataset datasets/expansion001 --until busco --dry-run
./run_dataset.sh submit --dataset datasets/expansion001 --until busco
./run_dataset.sh status --dataset datasets/expansion001
# After reviewing BUSCO and assembly outputs:
./run_dataset.sh submit --dataset datasets/expansion001 --until all
```

Endpoints are `assembly` (including CDS extraction), `busco`, `quant` (including
merge), `mapping`, and `all`. Earlier missing stages are included automatically.
Each upstream stage is a species array; only pending indices are submitted.
The species-to-index assignment remains fixed across stage arrays and retries.
The adapter translates it to GeneGalleon's sorted metadata filename index.
`array_size` defaults to 1000 and must be less than the site's `MaxArraySize`.
Species above that logical index are submitted in separate batches with a fixed
index offset, so a metadata table of thousands of species does not exceed the
scheduler's array index limit. Batches run in sequence; tasks within each batch
run in parallel up to `concurrency`. The submission receipt records the logical
indices and offset. Each task receives its own CPU, memory, and walltime allocation. Dependencies use `afterok`; failed
prerequisites cancel dependent jobs rather than publishing a partial dataset.

One dataset cannot be submitted again while its recorded jobs remain queued or
running. After failure, inspect the logs and rerun `submit`; completed species
products are reused. Status includes recorded worker state, job ID, and BUSCO
completeness. Logs and `submission_*.json` are under `datasets/<name>/jobs/`.
An interrupted/ambiguous scheduler submission is reported for reconciliation
before retry, to avoid submitting duplicate jobs. A timed-out assembler may have
to restart its unfinished stage; the adapter guarantees reuse of verified
completed stages, not checkpointing inside every external tool.

For a pilot without editing the master metadata:

```bash
./run_dataset.sh submit --dataset datasets/expansion001 --until quant \
  --species-list input/pilot_species.txt
```

Pilot submissions never trigger final dataset integration. Their completed
products are reused by the subsequent full submission. FASTQs and temporary work
are retained between stages. Array workers disable GeneGalleon's multispecies
summary; dataset-level collection is done once after all necessary results exist.
GeneGalleon provides its own locked shared reference preparation in its workspace;
all array tasks in a dataset share that workspace's downloads.

The downstream controller runs after the species arrays and uses the frozen
Slurm profile to schedule ODB chunks and other workflow jobs separately. Worker
resources come from the profile/rules, not the controller's CPU/memory allocation.
The controller also has a time limit: size it for the downstream queue/runtime
and resume the same dataset if necessary. A single species assembly exceeding
its own walltime still needs resource/time or assembly-setting adjustment.
No cluster submissions happen during tests or `--dry-run`.

## Removal and final integration

Delete species rows from the manual metadata and prepare a new dataset. The new
input view, expression, sequences, alignments, trees, and contrast pairs use only
that dataset's selected species. Previously completed species products and older
results are retained. Readding a species reuses the matching products.
Species-set-dependent analyses are recomputed for the new dataset; assembly,
BUSCO, quantification, and ODB annotations for unchanged species are reused.

Materialization requires all nonexcluded species to have validated CDS, BUSCO,
and quantification. An unsuccessful/missing species is not silently dropped.
BUSCO threshold exclusions are reported separately in `input_receipt.json`.
Prepared input views are published atomically and never scan all cached species.

`--until all` runs the configured downstream analyses and then the PhenoRadar
collector. Results are in `results/<name>/phenoradar_inputs/`. `--until mapping`
stops after the combined gene-to-OG mapping. Materialization can also be inspected
without launching downstream work:

```bash
./run_dataset.sh materialize --dataset datasets/expansion001
```

## Incremental ODB outside the dataset interface

The regular workflow also supports mixed mapping:

```yaml
odb:
  incremental: true
  existing_results: resources/odb_existing/tlight
  cache_dir: resources/odb_cache
  chunk_size: 20
  node: 3193
```

It validates selected protein hashes against explicit and cached snapshots,
imports matching species, maps only uncovered species, and publishes successful
new batches to the cache. A later dataset can reuse any subset of those batches.
OrthoDB version/node and known reference fingerprints must agree; changing the
reference requires a separate compatible snapshot/cache namespace. Legacy
snapshots with incomplete reference provenance retain that limitation. With
`incremental: false`, the original strict import/new-mapping behavior is retained.
