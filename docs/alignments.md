# OG protein alignments

[Back to README](../README.md)

This optional branch saves untrimmed amino-acid alignments for **every OG observed
in the selected species' ODB mappings**, retaining all mapped gene copies. It
reuses the translated protein FASTA and `orthogroups/mapping/mappings.sqlite`:

```text
Protein FASTA + ODB mappings -> collect by OG -> FAMSA -> untrimmed FASTA
```

## Running

To request alignments explicitly, run from the repository root:

```bash
./run_pipeline.sh --software-deployment-method conda \
  --configfile config/mydata.yaml --cores 16 --resources mem_gb=192 -- alignments
```

Prerequisite metadata selection, translation and ODB mapping are included
automatically. Existing mappings, including `odb.existing_results` snapshots,
are reused through the normal workflow rules. This target does not aggregate TPM,
annotate KO, or infer trees.

To include alignments in the default full workflow, set:

```yaml
alignment:
  enabled: true
  threads: 4
  mem_gb: 8
```

`enabled` defaults to `false`; the explicit `alignments` target works either way.
FAMSA 2.4.1 is supplied by `workflow/envs/alignment.yaml`. Threads are per OG;
`mem_gb` is the scheduling budget per branch job, including collection and
finalization. Concurrency follows the workflow's overall resource budget.

## Preserved information

- All copies from a species are kept, including identical sequences with different
  gene IDs. Multiple expression runs do not duplicate protein sequences.
- A gene assigned to several OGs appears in each. This is independent of
  `tpm.multimap`, which controls expression aggregation only.
- A singleton OG is saved directly without invoking FAMSA.
- OGs are not filtered by species count, copy count, sequence length, unknown
  residue fraction, or variation. No columns are trimmed or masked.
- Stop symbols (`*`, including terminal stops) and ambiguity symbols from the
  existing translation are preserved. No reading-frame correction or retranslation
  is performed here.

Inputs must contain nonempty, ungapped protein sequences (letters A–Z and `*`).
Mapped gene IDs must use `{species}_g{number}`, with an ASCII nonnegative integer
suffix and the exact metadata species ID. Malformed inputs, missing mapped genes,
gene IDs encoding a different species, and duplicate gene IDs fail the job
instead of silently dropping records. Multi-sequence outputs are checked for
unchanged gene IDs, equal aligned lengths and exact ungapped residue preservation.
These checks validate processing integrity, not orthology or alignment accuracy.

## Outputs and PhenoRadar

```text
results/<analysis>/orthogroups/alignments/
  {og}.faa
  provenance.json
```

FASTA headers retain the original gene IDs; output rows follow collected input
order. PhenoRadar's sequence inputs are just the `{og}.faa` files: the filename
gives the orthogroup, the first header token gives the gene ID, and removing the
final `_g{number}` gives the species. For example, `100007at3193.faa` containing
`>Abelia_chinensis_g0` identifies OG `100007at3193`, gene `Abelia_chinensis_g0`,
and species `Abelia_chinensis`. No `species=` attribute or `members.tsv` is added.
Species IDs match metadata and TPM outputs exactly, including underscores and
hyphens; only ODB protein filenames normalize hyphens. Splitting on the first
underscore is incorrect. Unmapped genes do not appear in these outputs.

PhenoRadar can use MSA column numbers directly. If needed, original protein
positions can be reconstructed by counting non-gap symbols along each row. No
site-coordinate table is stored. Column numbers refer to the particular saved
alignment and can change when its inputs or alignment settings change.

Representative-copy selection, OG/site selection and feature encoding are left
to PhenoRadar. A species without a sequence has no FASTA row; this alone is not
evidence of a biological gene deletion.

## Resuming and provenance

Collection reads each selected species' protein file once for sequence extraction,
keeps mappings for one species in memory, and caps simultaneously open OG files.
Collected FASTA and input records live in `work/<analysis>/orthogroups/alignments/inputs/`.
Each OG is a separate Snakemake job with a log and execution JSON under
`logs/<analysis>/orthogroups/alignments/` and a resource TSV in its `benchmarks/`
subdirectory. The execution JSON records the input/output hashes, actual
command, thread count and FAMSA executable hash; singletons record a direct copy.
The normal `run.json` records workflow code and resolved configuration.

After interruption, rerun the same command (with `--rerun-incomplete` if requested
by Snakemake). Completed OG jobs are reused. Changing abundance values or the TPM
ambiguity policy does not recompute alignments. Changes to selected species,
proteins or mappings rebuild collection and its dependent alignments.

Finalization verifies alignment hashes and publishes provenance only after all
current OG jobs have succeeded. It also removes a legacy `members.tsv` from its
output directory. On selection/mapping
updates it also removes obsolete `*.faa` files from the workflow-owned results
alignment directory, so the directory reflects the current OG set.
