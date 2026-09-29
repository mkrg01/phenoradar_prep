# OG protein alignments

[Documentation](index.md)

The `alignments` target uses FAMSA to align all mapped OGs, keeping every gene
copy without filtering or trimming.

## Running

```bash
./run_analysis.sh submit --analysis results/angiosperm_leaf_20260925/downstream/c4_photosynthesis_20260929 --target alignments
```

Set `alignment.enabled: true` before [preparing analysis](datasets.md#run-an-analysis)
to include alignments in `all`; the explicit target works either way. Each OG job
defaults to 4 CPUs/8 GB; see [resource overrides](running.md#resource-budgets).

## Input requirements

Gene IDs must be `{analysis_sample_id}_g{number}`: exact computational sample ID and nonnegative
integer. Proteins must be nonempty and ungapped (A–Z and `*`). Invalid inputs or
altered ungapped residues fail validation.

Identical copies, singletons, and stops are kept. Multi-OG genes appear in each
OG regardless of `tpm.multimap`. Each sample contributes its own assembled
proteins, including when several samples share a biological species.

## Outputs and PhenoRadar

Results are in `results/<build>/downstream/<analysis>/orthogroups/alignments/`:

- `{og}.faa`: untrimmed alignment with original gene IDs as FASTA headers.
- `provenance.json`: completed alignment inventory and provenance.

Remove the final `_g{number}` to recover the analysis sample ID, e.g.
`Abelia_chinensis_SRR14320411_g0` → `Abelia_chinensis_SRR14320411`. Use
`metadata/samples.tsv` to look up its biological species; do not split the sample ID.
Columns identify alignment sites; non-gap counts give protein positions.
PhenoRadar selects copies/sites. Missing sequences alone do not establish gene deletion.

Retries within the same analysis reuse completed OG jobs. New analysis IDs do
not share alignment results automatically. Logs: `results/<build>/logs/downstream/<analysis>/orthogroups/alignments/`.
