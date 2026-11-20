# OG protein alignments

[Documentation](index.md)

Align mapped OGs meeting a sample-occupancy threshold with FAMSA, retaining all
gene copies within selected OGs, then trim columns with trimAl in the same batch job.

## Running

```bash
./run_analysis.sh submit --analysis results/leaf/downstream/carnivory --target alignments
```

Set these conditions before preparing an analysis:

```yaml
alignment:
  enabled: true
  min_sample_fraction: 0.5
  trimal_mode: gappyout  # or automated1
```

Occupancy is the fraction of selected sample IDs with at least one mapped
protein in the OG. Separate runs of the same species count separately. A
threshold of `0.5` keeps groups present in at least half the samples. Use `0`
to align all observed groups or `1` to require every sample. The cutoff applies
only to OG alignments.

Each job uses four CPUs for one OG at a time. A batch holds up to 100 OGs and
requests 16 GB. An OG with an input FASTA of 16 MB or more gets its own 32 GB
job. Use `align_orthogroup_batch` and `align_large_orthogroup` to override
[resources](running.md#resource-budgets) for these jobs.

## Inputs and outputs

Proteins must be nonempty and ungapped, with gene IDs of the form
`{analysis_sample_id}_g{number}`. Identical copies, singletons, and stops are kept.
In this standalone target, multi-OG genes appear in each OG; expression
aggregation rejects such assignments.

Results are in `results/<build>/downstream/<analysis>/orthogroups/alignments/`:

- `{og}.faa.gz`: trimAl-filtered alignment with original gene IDs and all copies retained.
- `occupancy.tsv`: counts, fractions and selection decisions for all candidate OGs.
- `provenance.json`: completed alignment inventory and provenance.

The FASTA retains all mapped copies under their original gene IDs. trimAl
treats noncanonical residues as gaps while selecting columns and restores the
original characters in retained columns. Copies can become entirely gaps.
Positions refer to the trimmed alignment. For full-length coordinates, map a
reference protein back to that alignment.

A retry verifies and reuses matching alignments. Superseded alignments can
remain in the run's alignment cache and still occupy disk space. Use
[per-sequence TPM](phenoradar_inputs.md#amino-acid-inputs) to join expression by
sample and gene ID.
