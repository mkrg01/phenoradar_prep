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

Launcher tests exercise direct execution from another directory and simulate
Slurm's spooled script location and allocation settings. They check argument
forwarding, CPU and memory budgets, and exit-status propagation. A small real
Snakemake workflow verifies both modes without sbatch or srun submissions.
These tests do not require an actual Slurm allocation.

The automated tests do not validate real ODB assignments, reference download
availability, or Slurm execution. Run a real pilot with your input data to assess
mapping quality, runtime, and peak memory before launching the full analysis.
