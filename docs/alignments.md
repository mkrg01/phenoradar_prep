# OG protein alignments

[Documentation](index.md)

Align each mapped OG with FAMSA, retaining all gene copies without trimming.

## Running

```bash
./run_analysis.sh submit --analysis results/leaf/downstream/carnivory --target alignments
```

To include alignments in `all`, set `alignment.enabled: true` before preparing
the analysis. Weighted batches contain up to 500 OGs; inputs of 16 MB or more
run individually. Each batch defaults to 12 CPUs/32 GB, running up to three
FAMSA processes with four threads each; see
[resource overrides](running.md#resource-budgets).

OG alignments and provenance retain their original paths. A retried batch
verifies and reuses completed OGs; missing or changed members invalidate only
their batch receipt. `batches/*.json` records membership and reuse counts.
Batch scheduling uses input FASTA size as a cost estimate; exceptionally large
groups can still need a resource override.

## Inputs and outputs

Proteins must be nonempty and ungapped, with gene IDs of the form
`{analysis_sample_id}_g{number}`. Identical copies, singletons, and stops are kept.
In this standalone target, multi-OG genes appear in each OG; expression
aggregation rejects such assignments.

Results are in `results/<build>/downstream/<analysis>/orthogroups/alignments/`:

- `{og}.faa`: untrimmed alignment with original gene IDs.
- `provenance.json`: completed alignment inventory and provenance.

Remove the final `_g{number}` to recover the sample ID, then join
`metadata/samples.tsv` for biological species. Alignment columns identify sites;
non-gap counts give protein positions. PhenoRadar selects copies and sites.
Missing sequences alone do not establish gene deletion.
