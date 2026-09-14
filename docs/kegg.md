# KEGG annotation and KO expression

[Documentation](index.md)

The `kegg` target annotates translated proteins with KofamScan and sums original
input TPM by KEGG Orthology (KO), independently of ODB mapping.

## Run the branch

The branch prepares a missing `resources/kegg/snapshot_v1/` automatically.
See [KOfam/KEGG reference setup](references.md#kofam-and-kegg-reference) for
separate preparation, offline inputs, verification, and deliberate updates.

Add these settings to your dataset configuration:

```yaml
kegg:
  enabled: true
  threads: 4
  mem_gb: 8
  ambiguity: duplicate
```

`enabled: true` adds KEGG outputs to the default workflow. The explicit `kegg`
target works regardless of the flag:

```bash
./run_pipeline.sh --software-deployment-method conda \
  --configfile config/mydata.yaml --cores 24 --resources mem_gb=128 -- kegg
```

The rule-specific environment installs KofamScan 1.3.0 and its dependencies; its
executable is `exec_annotation`. Without Conda deployment, install KofamScan,
HMMER, GNU Parallel, and Ruby and expose them on `PATH`, along with the metadata
and translation dependencies.

Annotation runs once per species, searching all profiles in the snapshot.
Abundance-only or `kegg.ambiguity` changes repeat aggregation while reusing
annotations. Reuse checks bind results to proteins, references, code, executables,
and options. Failed searches retain diagnostics and retry fresh. Force the
annotation rule after replacing software in place, because Snakemake normally
schedules from file timestamps and recorded parameters.

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

Snakemake elapsed time and sampled resource use are saved separately at
`logs/<analysis>/kegg/benchmarks/<odb_species>.tsv`.

Per-run and support columns are `species`, `run`, `ko`, `tpm_sum`,
`annotated_genes`, and `quantified_genes`; the merged numeric long table contains
the first four columns. The counts refer to distinct genes assigned to each KO
under the selected policy and the subset represented in the run's abundance
table. Under `duplicate`, a multi-KO gene is counted once in each accepted KO.
Measured zeros remain zero. A KO with no quantified genes has blank `tpm_sum` in support
and wide tables and no row in the numeric long table. Partially quantified KOs
sum only observed genes, with counts exposing incomplete support. Merged tables
use the current manifest and omit stale species/runs from earlier selections.

Support and MODULE membership do not establish genomic absence, pathway activity,
flux, or reaction/MODULE completeness.

## Passing KO features to PhenoRadar

Use the [input collector](phenoradar_inputs.md) to link KO expression and the
requested KO-to-module/pathway maps. It requires one run per species and checks
expression coverage. The collection guide gives the KO column settings.

## References

- [KofamScan documentation](https://github.com/takaram/kofam_scan)
- [KofamKOALA paper](https://doi.org/10.1093/bioinformatics/btz859)
- [KEGG API manual](https://www.kegg.jp/kegg/rest/keggapi.html)
- [KEGG MODULE definitions](https://www.kegg.jp/kegg/module.html)
