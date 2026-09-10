# KEGG annotation and KO expression

[Back to README](../README.md)

This optional branch annotates translated proteins with KofamScan and sums the
original abundance-table TPM by KEGG Orthology (KO). It runs independently of ODB
mapping. Preserved gene IDs connect protein annotations to `target_id` values;
OG TPM is never used to reconstruct KO expression.

## Prepare a reference snapshot

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

## Run the branch

Add these settings to your dataset configuration:

```yaml
kegg:
  enabled: true
  threads: 4
  mem_gb: 8
  ambiguity: duplicate
```

`enabled: true` adds KEGG outputs to the default full workflow. With the default
`false`, OG-only runs require no KOfam reference or software. The explicit
`kegg` target requests this branch regardless of the flag, without ODB mapping:

```bash
./run_pipeline.sh --software-deployment-method conda \
  --configfile config/mydata.yaml --cores 24 --resources mem_gb=128 -- kegg
```

The rule-specific environment installs KofamScan 1.3.0 and its dependencies; its
executable is `exec_annotation`. Without Conda deployment, install KofamScan,
HMMER, GNU Parallel, and Ruby and expose them on `PATH`. Standard metadata and
translation dependencies are also required. The resource values above are
starting budgets, not measured requirements for every dataset.

Annotation runs once per species and is reused across its runs. Abundance-only
updates repeat TPM aggregation without repeating annotation. When an annotation
job runs, its content fingerprint binds reuse to the protein, reference, code,
executable, and options; completed outputs are verified before reuse. Failed
searches retain diagnostic work and retry fresh. Partial output is never treated
as complete. Snakemake normally schedules jobs from file timestamps and recorded
parameters; force the annotation rule after replacing software in place without
changing its configured path.

Changing only `kegg.ambiguity` also repeats aggregation and merging while reusing
completed species annotations. Explicit `drop` settings in older dataset
configurations continue to exclude multi-KO genes; omit the setting or choose
`duplicate` to use the new default.

All profiles in the snapshot are searched. For deliberately restricted sets,
record the choice and assess coverage, including organellar genes, before
interpreting unassigned KOs biologically.

## Assignment and quantification

For HMMER input, one terminal stop from CDS translation is removed and recorded
per gene. Internal stops are rejected. Original IDs and source FASTA files are
preserved. Raw `detail.tsv` retains scores, thresholds, and E-values. Acceptance
uses KofamScan's marker, calculated before printed scores are rounded. Unmarked
hits and hits without a defined threshold are not quantitative assignments.
An HMM score is not a probability of correct annotation.

| Gene status | Meaning |
| --- | --- |
| `unique` | One distinct accepted KO; eligible for TPM aggregation |
| `ambiguous` | Multiple accepted KOs; full TPM contributes to each by default |
| `threshold_missing` | No accepted KO and at least one hit without a threshold |
| `below_threshold` | Hits exist but none passes its threshold |
| `unannotated` | No reported KO hit |

`kegg.ambiguity: duplicate` is the default: each gene contributes its full
original TPM to every distinct accepted KO. For example, a gene with TPM 10
and two accepted KOs adds 10 to each KO, for a sum of 20 across those KOs.
Repeated hits for the same gene/KO pair do not add extra contributions.
Below-threshold and threshold-missing hits never receive TPM.

`drop` excludes multi-KO genes from KO sums. `error` rejects a run with a
quantified multi-KO gene, including an observed zero. Candidate annotations
are retained under every policy. The annotation status `ambiguous` records
multiple accepted KOs; it does not distinguish true multifunctionality from
alternative annotation candidates. Equal splitting is not implemented.
In `genes.tsv`, `selected_ko` is populated only for a unique assignment.
Use the `accepted=1` rows in `gene_kos.tsv` to retrieve all accepted KOs,
including those contributing under `duplicate`.

KO values sum **original input TPM without renormalizing** to the retained KO
set. Under `duplicate`, the sum across KOs can exceed the total input TPM;
these KO features overlap and their sum is not a transcript total.
QC reports total TPM, retained TPM/fraction, TPM by annotation status,
and targets without translated proteins. `retained_targets` and `retained_tpm`
count each contributing input gene once, so `retained_tpm_fraction` remains
between zero and one. `quantified_assignments` counts contributing gene/KO pairs,
and `ko_tpm_sum` sums all KO contributions, including repeated TPM across KOs.
`ambiguous_tpm` counts the original TPM of multi-KO genes once; it is included
in retained TPM under `duplicate` and excluded under `drop`.
Runs with no retained KOs still produce QC and a header-only numeric expression
table. Nonempty all-zero
abundance tables are valid in this branch; zero-denominator fractions are blank
in TSV and null in JSON. Duplicate target IDs and negative/nonfinite values are
errors.

## Outputs

```text
results/<analysis>/kegg/
  reference_qc.json
  ko_modules.tsv                   # ko, module
  ko_pathways.tsv                  # ko, pathway
  species/<odb_species>/
    detail.tsv                     # Raw KofamScan output
    gene_kos.tsv                   # Candidate hits, thresholds, acceptance
    genes.tsv                      # Every protein, including unannotated genes
    provenance.json
    benchmark.tsv                  # Snakemake elapsed seconds and sampled resource use
    execution_config.json
    stdout.log
    stderr.log
  runs/<run>.tsv                   # Per-run expression/support
  runs/<run>.qc.json
  gene_kos.tsv                     # Hits for selected species
  genes.tsv                        # Genes for selected species
  ko_tpm_sum.tsv                   # Numeric long expression
  ko_tpm_sum_wide.tsv              # One row per species/run
  ko_support.tsv                   # All annotated KO/run support records
  mapping_qc.tsv                   # One row per selected run
```

Per-run and support columns are `species`, `run`, `ko`, `tpm_sum`,
`annotated_genes`, and `quantified_genes`; the merged numeric long table contains
the first four columns. The counts refer to distinct genes assigned to each KO
under the selected policy and the subset represented in the run's abundance
table. Under `duplicate`, a multi-KO gene is counted once in each accepted KO.
Measured zeros remain zero. A KO with no quantified genes has blank `tpm_sum` in support
and wide tables and no row in the numeric long table. Partially quantified KOs
sum only observed genes, with counts exposing incomplete support. Merged tables
use the current manifest and omit stale species/runs from earlier selections.

Support is not proof of genomic absence or reaction completeness. MODULE
membership does not capture required steps and alternative reactions. This
extension does not compute metabolic flux, pathway activity, or MODULE
completeness.

## Passing KO features to PhenoRadar

After selecting compatible tissue/condition runs and explicitly constructing a
species-level table, KO identifiers can be used as feature IDs:

```yaml
data:
  tpm_path: ko_expression_species.tsv
  feature_col: ko
  value_col: tpm_sum
  orthogroup_annotation_path: null
```

Prep preserves runs; it does not choose or average them into a species. Do not
pass multiple runs of one species to a consumer that sums duplicate
species/feature rows. Check support/QC before using consumers whose default
missing-feature value is zero.

Module scoring, learned scaling, phenotype-dependent selection, and group-lasso
fitting belong inside PhenoRadar's training folds. These outputs provide fixed,
label-independent annotation and expression inputs.

## References

- [KofamScan documentation](https://github.com/takaram/kofam_scan)
- [KofamKOALA paper](https://doi.org/10.1093/bioinformatics/btz859)
- [KEGG API manual](https://www.kegg.jp/kegg/rest/keggapi.html)
- [KEGG MODULE definitions](https://www.kegg.jp/kegg/module.html)
