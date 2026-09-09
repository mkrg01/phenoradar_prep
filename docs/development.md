# Testing and validation

[Back to README](../README.md)

Run the tests from the repository root. For workflow execution, see the [running guide](running.md).

Tests require a Python environment containing `pytest`, `pandas`, `ete4`,
`matplotlib`, and `PyYAML`. The integration test also requires Snakemake and seqkit.
Provide their executable paths if they are not on `PATH`:

```bash
SNAKEMAKE_BIN=/path/to/snakemake \
SEQKIT_BIN=/path/to/seqkit \
  /path/to/test-environment/bin/python -m pytest -q tests
```

The suite covers selection, multiple runs, invalid inputs, taxonomy snapshots,
ODB output formats, checksum verification, duplicate and ambiguous mappings, TPM
normalization, and recovery after failure. The integration test uses real seqkit
and a test-only ODB substitute. It checks a complete workflow, separate preparation
followed by execution, unchanged reruns, abundance updates, and reduced species
selections. Without Snakemake or seqkit, that integration test is skipped.

Taxonomy bootstrap tests build a real ETE4 database from a small synthetic
taxdump, substituting only the network response. They check download/build
failure recovery, local source copying, and reuse of existing snapshots. The
workflow integration test starts with a missing taxonomy snapshot and verifies
that later runs no longer require the bootstrap source. Tests do not download
the full NCBI taxonomy.

Launcher tests exercise direct execution from another directory and simulate
Slurm's spooled script location and allocation settings. They check argument
forwarding, CPU and memory budgets, and exit-status propagation. A small real
Snakemake workflow verifies both modes without sbatch or srun submissions.
These tests do not require an actual Slurm allocation.

KEGG tests additionally cover frozen reference snapshots and corruption checks,
KofamScan detail-TSV interpretation, ambiguity, original-TPM accounting, missing
observations, observed zeros, and runs with no retained KOs. A second workflow
integration test uses a test-only KofamScan substitute with real Snakemake and
seqkit. It checks the standalone `kegg` target, opt-in full workflow, annotation
reuse across runs, abundance-only updates, and removal of stale species from
merged outputs. No reference downloads are needed for these tests.

The automated tests do not validate real ODB assignments, reference download
availability, or Slurm execution. Run a real pilot with your input data to assess
mapping quality, runtime, and peak memory before launching the full analysis.
Likewise, synthetic KofamScan tests do not establish biological KO assignment
accuracy or calibrate confidence for novel plant sequences.
