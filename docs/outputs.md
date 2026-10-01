# Outputs and QC

[Documentation](index.md) · [PhenoRadar inputs](phenoradar_inputs.md)

## Directory layout

| Path | Contents |
| --- | --- |
| `results/<build>/database/` | Portable sample products and provenance |
| `results/<build>/downstream/<analysis>/` | Analysis settings, inputs, and results |
| `results/<build>/work/` | GeneGalleon workspaces, computation, and build caches |
| `results/<build>/logs/` | Rule logs and benchmarks for the database and analyses |
| `results/<build>/jobs/`, `<analysis>/jobs/` | Submission records, stage status, and launcher logs |
| `resources/` | Shared reference and software caches |

Keep generated products unchanged: published databases may share hard links with
working files. To move a database, [copy the complete bundle](datasets.md#copying-a-completed-build-to-another-project).

## Database files

Paths below are relative to `results/<build>/database/`. Filename placeholders
follow the [sample identifier conventions](inputs.md#identifiers).

| Path | Contents |
| --- | --- |
| `metadata.tsv` | Metadata after run exclusions |
| `cds/{species}_longestCDS.fa.gz` | Longest CDS per gene |
| `busco/summary.tsv`, `busco/full/` | Per-sample counts and full BUSCO tables |
| `quant/{species}/{run}/{run}_abundance.tsv` | Gene expression: `target_id`, `tpm` |
| `proteins/{odb_species}_protein.fa` | Translated proteins |
| `odb/snapshot.json`, `odb/species/*.tsv.gz` | Mapping provenance and gene/OG tables |
| `expression/runs/{run}.tsv`, `{run}.qc.json` | Saved OG expression and QC |

## Result files

```text
results/<build>/downstream/<analysis>/
  metadata/                  # Selected samples, traits, and BUSCO QC
  orthogroups/
    mapping/                 # gene_id / orthogroup tables and provenance
    expression/
      tpm_sum.tsv            # Original TPM summed by OG
      tpm.tsv                # Rescaled OG TPM
      tpm_sum_wide.tsv
      tpm_wide.tsv
      mapping_qc.tsv
    alignments/              # Optional all-copy alignments
  kegg/                      # Optional KO annotation and expression
  phylogeny/                 # Trees, dating, pairs, and taxonomy reports
  phenoradar_inputs/         # Collected inputs for PhenoRadar
```

Start with these checks:

| File | Review |
| --- | --- |
| `metadata/selection.json` | Which samples were included and why |
| `metadata/samples.tsv` | Sample IDs, biological species, and runs |
| `metadata/busco_completeness.svg` | Per-sample BUSCO completeness |
| `orthogroups/expression/mapping_qc.tsv` | Mapped/retained TPM fractions and ambiguous targets |

Optional [analysis guides](index.md#choose-an-analysis) describe their own outputs.
Keep `run.json`, branch provenance, reference snapshots, and the release's
`image.json` with your results.

## Storage cleanup

ODB, KofamScan, and GeneGalleon remove temporary files after saving and verifying
their results. ODB cleans up each completed chunk. Downloaded FASTQ/SRA files are
removed after both assembly and quantification finish. Original inputs, reusable
results, and failed work are kept. For debugging, [retain intermediates](configuration.md#intermediate-storage)
when preparing the run.

Check cleanup status and preview files left by completed jobs:

```bash
./run_build.sh status --build results/leaf --storage
./run_build.sh cleanup --build results/leaf
./run_analysis.sh cleanup --analysis results/leaf/downstream/carnivory
```

Add `--apply` to a `cleanup` command to delete the listed files. Cleanup errors
are reported as `pending`; retrying cleanup does not rerun computation.
Keep the rest of `work/` for reuse.

## TPM interpretation

Long tables contain `species` (sample ID), `run`, `orthogroup`, and `tpm_sum` or
`tpm`. Wide tables have one row per run and one column per OG; missing combinations
are zero.

- **`tpm_sum`** sums original gene TPM within each OG, excluding unmapped genes.
- **`tpm`** rescales retained OG values to one million per run. It describes
  relative expression within the retained OG set.

Duplicate gene/OG pairs count once. A quantified gene assigned to multiple OGs
stops aggregation. Duplicate targets, negative/nonfinite TPM, and no positive
retained TPM also cause errors.

[PhenoRadar export](phenoradar_inputs.md) removes the `run` column and preserves
TPM values; it never averages samples. [KO expression](kegg.md#outputs) uses
original TPM sums and distinguishes unquantified features from measured zeros.
