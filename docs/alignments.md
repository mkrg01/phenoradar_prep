# OG protein alignments

[Documentation](index.md)

Align each mapped OG with FAMSA, retaining all gene copies without trimming.

## Running

```bash
./run_analysis.sh submit --analysis results/leaf/downstream/carnivory --target alignments
```

To include alignments in `all`, set `alignment.enabled: true` before preparing
the analysis. Each OG job defaults to 4 CPUs/8 GB; see
[resource overrides](running.md#resource-budgets).

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
