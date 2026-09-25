# OG protein alignments

[Documentation](index.md)

The `alignments` target uses FAMSA to align all mapped OGs, keeping every gene
copy without filtering or trimming.

## Running

```bash
./run_pipeline.sh --cores 16 --resources mem_gb=192 --configfile analyses/analysis001/pipeline.yaml -- alignments
```

Prerequisites run automatically. Set `alignment.enabled: true` to include this
in `all`; see [resources](running.md#resource-budgets) for job sizing.

## Input requirements

Gene IDs must be `{species}_g{number}`: exact metadata species ID and nonnegative
integer. Proteins must be nonempty and ungapped (A–Z and `*`). Invalid inputs or
altered ungapped residues fail validation.

Identical copies, singletons, and stops are kept. Multi-OG genes appear in each
OG regardless of `tpm.multimap`; expression replicates do not duplicate sequences.

## Outputs and PhenoRadar

Results are in `results/<run_name>/orthogroups/alignments/`:

- `{og}.faa`: untrimmed alignment with original gene IDs as FASTA headers.
- `provenance.json`: completed alignment inventory and provenance.

Remove the final `_g{number}` to recover the species, e.g. `Abelia_chinensis_g0`
→ `Abelia_chinensis`. Columns identify alignment sites; non-gap counts give
protein positions. PhenoRadar selects copies/sites. Missing sequences alone do
not establish gene deletion.

Rerunning reuses completed OG jobs. Abundance/TPM-policy changes leave alignments
intact; sequence/membership changes rebuild them, potentially changing columns
and removing obsolete OG files. Logs: `logs/<run_name>/orthogroups/alignments/`.
