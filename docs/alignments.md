# OG protein alignments

[Documentation](index.md)

Align mapped OGs meeting a sample-occupancy threshold with FAMSA, retaining all
gene copies within selected OGs without trimming.

## Running

```bash
./run_analysis.sh submit --analysis results/leaf/downstream/carnivory --target alignments
```

Set these conditions before preparing an analysis:

```yaml
alignment:
  enabled: true
  min_sample_fraction: 0.5
```

The default `0.5` selects OGs represented in at least 50% of the selected sample
IDs in `metadata/samples.tsv`, after sample exclusions and BUSCO selection.
The required count is rounded up: 5,582 samples require at least 2,791.
Each sample counts once even with multiple gene copies; distinct sample IDs from
the same biological species count separately. All copies of a retained OG are
aligned. Allowed values are 0–1; `0` retains all observed OGs and `1` requires
every selected sample. An empty qualifying set produces an empty inventory.
This condition affects OG alignments, not TPM aggregation or BUSCO-marker trees.

Weighted batches contain up to 500 selected OGs; inputs of 16 MB or more
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
- `occupancy.tsv`: counts, fractions and selection decisions for all candidate OGs.
- `provenance.json`: completed alignment inventory and provenance.

After successful completion, the directory's `*.faa` files are exactly the
selected OGs listed in `provenance.json`. Former outputs below the threshold are
moved to `work/<run>/orthogroups/alignments/cache/` (under `work_root` when set).
Re-adopted OGs are verified against their input, alignment and aligner records,
then moved back without recomputation. Deleting this optional cache saves disk
space but requires recomputation if those OGs are selected again.
OGs no longer present in the current raw collection are removed from both the
final directory and cache, as with obsolete outputs before this cache existed.

Moving excluded results shrinks the final output directory; retaining the cache
does not reduce total stored alignment data. Files are moved rather than kept
in both locations. Across filesystems, a temporary copy is needed during a move.

Remove the final `_g{number}` to recover the sample ID, then join
`metadata/samples.tsv` for biological species. Alignment columns identify sites;
non-gap counts give protein positions. PhenoRadar selects copies and sites.
Missing sequences alone do not establish gene deletion.
