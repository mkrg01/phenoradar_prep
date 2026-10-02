# Angiosperm leaf

One leaf RNA-seq run per original NCBI taxid, searched with:

```text
"Angiosperms"[Organism] AND "Illumina"[Platform] AND "RNA-seq"[Strategy]
```

See the [shared layout and selection policy](../README.md).

## Initial state

This directory initially contains settings, rules, and historical exclusions.
`metadata.tsv`, `accepted_samples.tsv`, `selection.tsv`, and `provenance.json`
are created only when the workflow is run. Begin with a new `amalgkit metadata`
query covering both historical and newly registered runs; then curate that new
download and choose representatives. Raw metadata and candidates go under
`work/`; `update` creates the dataset's `metadata.tsv`. Use `--dry-run` for a
candidate-only preview when desired.
Missing metadata and success tables mean there is no previous selection to
retain, so no initialization command is required.
The 744 historical exclusions are enriched with their original `taxid` and
`bioproject`. Their reasons and source build are preserved. Historical metadata
was consulted only to resolve those exclusion identities during migration.

The old database was produced with older software outside the fixed container.
Consequently, no success evidence is provided. `build.yaml` sets `reuse_from:
auto` to reuse compatible, verified sample stages from builds under `results/`.
The old database is not automatically imported. Newly selected runs without
matching stage receipts are rebuilt. After QC review, use `record` to adopt
successful representatives, including those from a partially completed build.

The GeneGalleon image is fixed by OCI digest, version 0.7.77, and source commit
in `build.yaml`. Its bundled AMALGKIT is **0.16.76**, commit
`e004851c63f0f902aa48f69e270ba42bf2d2f1ba`. It is the pinned version used for
preprocessing, rather than an independently installed AMALGKIT version.

The effective rules start from this image's 711 plantae rules, with the five
historical additions and leaf-only parameter carried forward: 716 rules total.
The additions recognize fronds as leaves and exclude single-cell, spatial,
cell-type/LCM, and special RNA assays. Their origins and migration are recorded
in `rules/legacy_origin.json` and `rules/initial_migration.json`; the old
effective table is retained under `rules/legacy/` for comparison.

## Updating metadata

Choose between a fresh NCBI query, reselection from a saved curated source, and
processing retries using the [restart entry points](../README.md#choosing-where-to-restart).
Metadata updates use a new attempt directory and a new build name. Processing
retries use the existing build and keep its selected runs unchanged.

The host Python needs PyYAML (use the provided metadata environment).
Use `WORKFLOW_PYTHON` to select another interpreter if needed. Apptainer or
Singularity runs the pinned image; AMALGKIT does not need a host installation.

Run from the repository root. Prepare the host environment and pinned image:

```bash
conda env create -n phenoradar_metadata -f datasets/environment.yaml
conda activate phenoradar_metadata
export WORKFLOW_PYTHON="$(command -v python)"
./run_build.sh fetch-software --config datasets/angiosperm_leaf/build.yaml
```

`fetch-software` downloads the pinned GeneGalleon source and SIF into the software
cache. The subsequent metadata commands use this image.

Edit `rules/select_rules.tsv` and `excluded_accessions.tsv` before curation.
The upstream rule snapshot is supplied with the dataset definition. After deliberately
changing container pins, refresh it with:

```bash
./run_metadata.sh rules --refresh-upstream
```

This replaces only `rules/upstream/`. Compare it against the effective rules and
merge any desired changes yourself; your edited rules are never overwritten.

Choose a new attempt directory, then update metadata in one command:

```bash
attempt=work/datasets/angiosperm_leaf/20261001
./run_metadata.sh update --work "$attempt"
```

`update` performs the NCBI query, runs AMALGKIT select with saved rules, chooses
representatives, and writes the dataset metadata. Its fetch and curate stages
save the image checksum, AMALGKIT
version/commit, command, input hashes, and `amalgkit.log`. Each stage refuses to
overwrite an existing attempt. These operations do not submit analysis jobs.
Use `update --work "$attempt" --dry-run` instead to stop at a candidate preview.
The individual `fetch`, `curate`, `select`, and `accept` commands remain available.

After changing rules, reuse the **newly fetched** raw table to retry curation
without querying NCBI again. Choose a new output attempt:

```bash
retry=work/datasets/angiosperm_leaf/20261001_rules_v2
./run_metadata.sh curate --work "$retry" \
  --metadata "$attempt/fetch/metadata/metadata.tsv"
./run_metadata.sh update --work "$retry" \
  --metadata "$retry/curate/metadata/metadata.tsv"
```

An update normally starts with a new `fetch`. Rule adjustments within that
update can use its frozen download. After changing only exclusions or overrides,
use `update --work <new-attempt> --metadata <previous-curated.tsv>` to select again
without another NCBI query. `select` records the curated source's
checksum and provenance. In this pinned AMALGKIT version, `sample_group` is
configured in `select_rules.tsv`; `amalgkit select --sample_group leaf` is not a
supported CLI option.

Inspect `candidate/metadata.tsv`, `candidate/selection.tsv`, and
`candidate/provenance.json`. Selection reports distinguish retained successes,
reselected historical runs, replacements, new taxids, and taxids with no eligible
run, with previous/selected BioProjects and `total_bases` for comparison.
To change a choice, edit effective rules or add a TSV row to `overrides.tsv`:

```tsv
taxid	run	reason
12345	SRR12345678	Verified leaf tissue in the original study
```

An override must reference an eligible candidate; it cannot bypass an exclusion.
Regenerate into a new candidate directory after editing. Candidate metadata and
reports are checked against their hashes on acceptance, so use overrides rather
than editing the generated table directly. Rules, exclusions, successful-run
evidence, overrides, and build settings must also match those used for selection.

If you used `--dry-run`, publish the inspected candidate with:

```bash
./run_metadata.sh accept --candidate "$attempt/candidate"
```

This updates the dataset's metadata, selection report, and provenance. New runs
still have no success evidence. Retained runs with source changes have
informational review flags; they remain selected until explicitly excluded or
overridden. Exclusions match run accessions only, so another eligible run with
the same taxid and BioProject can replace an excluded representative.

## Running the database build separately

Activate the [normal workflow environment](../../docs/running.md#installation-and-normal-execution)
and set `WORKFLOW_PYTHON` to its Python interpreter before running build commands.
Use the dataset's build config explicitly; the repository's default build config
has not been changed. No analysis starts during metadata preparation.

```bash
./run_build.sh plan --config datasets/angiosperm_leaf/build.yaml
./run_build.sh prepare --config datasets/angiosperm_leaf/build.yaml \
  --name angiosperm_leaf_20261001
```

`plan` shows required work. `prepare` saves reviewed inputs and settings without
submitting jobs. Start processing when ready with:

```bash
./run_build.sh submit --build results/angiosperm_leaf_20261001
```

Results go under `results/angiosperm_leaf_20261001/`; the portable database goes
under its `database/` directory. Nothing is published automatically.

## Recording new successes for the next update

Inspect completed sample stages, their BUSCO counts, and proposed exclusions,
even while other samples are incomplete:

```bash
./run_metadata.sh record --build results/angiosperm_leaf_20261001 --dry-run
```

After reviewing QC, add unusable runs to `excluded_accessions.tsv`. An accession
and optional reason suffice; taxid and BioProject are optional context. Then
record QC for current representatives:

```bash
./run_metadata.sh record --build results/angiosperm_leaf_20261001
# To record QC only for a reviewed subset:
./run_metadata.sh record --build results/angiosperm_leaf_20261001 --runs SRR123 SRR456
```

When adding rows to the supplied exclusion table, leave unused context cells
blank while keeping the table's existing columns.

`record` requires matching pinned-container provenance and verified products.
It skips already excluded and superseded runs. New samples with completed
assembly and BUSCO products below `busco_threshold` are reported as
`busco_below_threshold` and automatically added to `excluded_accessions.tsv`,
even if quantification is not finished. The recorded reason is
`busco_completeness_below_<threshold>`, with taxid, BioProject, and source build.
Existing rows, manual reasons, and custom columns are preserved; missing context
columns are added. Repeating the command does not duplicate exclusions.

New adoption requires completed assembly, BUSCO, and quantification products and
`(BUSCO single + duplicated) / total >= busco_threshold`. `selection.yaml` sets
the threshold to `0.5`; exactly 50% passes. Existing adoption records are preserved
unless explicitly excluded; later threshold changes do not revoke them or add
automatic exclusions for those adopted runs. Missing or failed assembly/BUSCO
stages are reported as incomplete and remain available for retry.

`record` updates the success and exclusion tables but leaves metadata unchanged.
`--dry-run` only previews, and `--runs` limits both adoption and automatic
exclusion to the reviewed subset. Invocation without `--dry-run` records these QC
decisions. `baseline` imports a reviewed completed database with the same
adoption threshold and records its below-threshold runs as exclusions.

After recording QC, select replacement runs using the saved curated source and
a new attempt directory:

```bash
./run_metadata.sh update --work work/datasets/angiosperm_leaf/reselection_01 \
  --metadata "$attempt/curate/metadata/metadata.tsv"
```

For each excluded run's taxid, this selects the eligible run with the greatest
`total_bases`, including candidates in the same BioProject. A taxid with no
eligible candidate leaves the metadata. Prepare and submit a new build name to
process the replacements; `reuse_from: auto` reuses other completed sample stages.

For a sample defect, edit the exclusions, run `update` with a saved curated table,
and prepare a new build name. With `reuse_from: auto`, successful sample stages
are reused from previous builds, including unfinished ones. Prepared builds
keep their original inputs. For memory or time failures, leave exclusions alone
and resubmit the same build with `--resources <retry.yaml>`; successful stages
remain reusable. Automatic stage reuse does not reuse unfinished assembly
scratch or mapping/expression caches. Explicit completed-database reuse can
also reuse mapping and expression products.
