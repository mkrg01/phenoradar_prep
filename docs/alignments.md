# OG protein alignments

[Documentation](index.md)

The `alignments` target uses FAMSA to create untrimmed protein alignments for
all mapped OGs, retaining every gene copy.

## Running

```bash
./run_pipeline.sh --cores 16 --resources mem_gb=192 -- alignments
```

Selection, translation, and ODB mapping are scheduled as needed. To include this
branch in `all`:

```yaml
alignment:
  enabled: true
```

Each `align_orthogroup` job aligns one OG across species, retaining all copies,
with a default of 4 threads and 8 GB total memory. Sequence collection and
finalization each use 1 thread and 8 GB. See
[resource budgets](running.md#resource-budgets) for concurrency and per-rule overrides.
The explicit target also works when disabled.

## Preserved information

All copies are kept, including identical sequences with different IDs. A gene
mapped to several OGs appears in each, independently of `tpm.multimap`. Replicate
expression runs do not duplicate sequences. Singletons are copied directly.
There is no OG/site filtering or trimming; translation symbols, including stops,
are preserved.

Mapped gene IDs must be `{species}_g{number}` using the exact metadata species ID
and a nonnegative integer. Inputs must be nonempty, ungapped proteins (A–Z and `*`).
Missing genes, invalid IDs, and changed ungapped residues fail validation.

## Outputs and PhenoRadar

```text
results/<run_name>/orthogroups/alignments/
  {og}.faa
  provenance.json
```

The filename identifies the OG; FASTA headers preserve gene IDs. Remove the final
`_g{number}` to recover the species, including its underscores and hyphens.
For example, `>Abelia_chinensis_g0` belongs to `Abelia_chinensis`.

Alignment columns identify sites; count non-gap residues for protein positions.
Columns can change after input/settings updates. PhenoRadar handles copy and
site selection. A missing sequence alone does not establish gene deletion.

## Resuming and provenance

Rerun the same command to reuse completed OG jobs. Abundance-only and TPM-policy
changes do not recompute alignments. Species/protein/mapping changes rebuild them;
obsolete OG FASTAs are removed from the workflow-owned results directory.

Logs and per-OG execution records are under
`logs/<run_name>/orthogroups/alignments/`; `provenance.json` records the completed set.
