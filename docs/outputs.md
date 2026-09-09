# Outputs and TPM interpretation

[Back to README](../README.md)

Output paths use the `analysis` name from your [configuration](configuration.md).

## Output files

```text
results/<analysis>/
  run.json                          # Resolved configuration, code hashes, Python environment
  metadata/
    metadata_all.tsv
    metadata_high_busco.tsv
    samples.tsv
    species_high_busco.txt
    selection.json                  # Selection counts and input/taxonomy hashes
    busco_completeness.svg
  proteins/                         # Species FASTA files and translation provenance
  odb/
    manifests/                      # Chunk plan and FASTA manifests
    chunks/chunk_000/                # Annotations, hits, summary, provenance, native results
    merged/
      gene_orthogroups.tsv           # Unique #query / ODB_OG pairs
      mappings.sqlite               # Indexed gene-to-OG mappings
      merge_qc.json
  tpm/
    runs/                           # Per-run results and QC JSON
    tpm_sum.tsv                     # Sums of input TPM by orthogroup
    tpm.tsv                         # OG TPM rescaled to one million per run
    tpm_sum_wide.tsv
    tpm_wide.tsv
    mapping_qc.tsv
```

Logs and ODB resource benchmarks are saved under `logs/<analysis>/`. Temporary
ODB work is stored under `work/<analysis>/odb/`.

## TPM interpretation

Long tables contain `species`, `run`, `orthogroup`, and either `tpm_sum` or `tpm`.
Wide tables contain one row per run, identified by `species` and `run`, with one
column per orthogroup. Missing run/orthogroup combinations are filled with zero.
Multiple runs from a species are not pooled or averaged.

- `tpm_sum` sums the original input TPM values assigned to each orthogroup,
  without rescaling. Unmapped genes are excluded.
- `tpm` rescales retained orthogroup values to sum to one million within each run.
  It therefore describes relative expression within the retained OG set.

Duplicate gene/OG pairs are removed before aggregation. If a gene maps to multiple
orthogroups, `tpm.multimap` controls its treatment:

| Policy | Behavior |
| --- | --- |
| `error` (default) | Stop and report the ambiguous genes |
| `drop` | Exclude genes assigned to multiple OGs |
| `split` | Divide each gene's TPM equally among its assigned OGs |

`mapping_qc.tsv` reports the fraction of input TPM mapped to OGs, ambiguous target
counts, the retained TPM fraction, and other mapping statistics. Aggregation
rejects duplicate target IDs, negative or nonfinite TPM values, and runs with no
positive TPM retained after mapping and ambiguity handling.

## Optional KEGG outputs

The [KEGG branch](kegg.md) writes `results/<analysis>/kegg/`, separately from the
OG tables above. Its `ko_tpm_sum.tsv` contains sums of original input TPM and is
**not renormalized** to the retained KO set. `ko_support.tsv` records annotated
and quantified gene counts, including KOs with no quantified genes. Such KOs have
blank values in support/wide tables and are omitted from the numeric long table;
observed zero expression remains zero. Run identities are preserved.
