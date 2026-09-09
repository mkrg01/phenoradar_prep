# Testing and validation

[Back to README](../README.md)

Run the tests from the repository root. For workflow execution, see the [running guide](running.md).

Tests require a Python environment containing `pytest`, `pandas`, `ete4`,
`matplotlib`, `PyYAML`, `numpy`, and `biopython`. The integration test also requires Snakemake and seqkit.
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

[Phylogeny tests](phylogeny.md#outputs-and-validation) additionally exercise
BUSCO identifier mapping, selection by highest overall coverage with deterministic
ID ties, independence from order labels, replicate counts and BUSCO match lengths,
retention of low-coverage loci before and after alignment QC,
cdskit padding/masking/translation, alignment QC with low/zero informative-site
counts, retention of species present in only one gene tree, missing-species
rejection and count diagnostics, and calibration consistency. Install the pinned cdskit
from `workflow/envs/phylogeny.yaml` in the test environment; tests compare its
actual CLI with the workflow's Python API on stops, partial codons, frame changes,
IUPAC/gap codons and alternate genetic codes. When tool paths are supplied, they
run real FAMSA, trimAl (`TRIMAL_BIN`), VeryFastTree, ASTRAL-IV and treePL on
synthetic inputs. They verify selected columns, original-X retention, unchanged
reruns and FAMSA reuse after changing trimAl mode. They check that dating
uses time units, honors calibrations, and reuses the existing species tree.
The treePL tests also exercise repeated random CV on a 32-tip synthetic tree
with heterogeneous lengths and a zero branch, reject incomplete/nonfinite
outputs, and check rooted topology and ultrametricity. The project build includes
the documented [upstream fixes](../workflow/patches/README.md); set `TREEPL_BIN`
to that executable. A separate synthetic smoke run additionally exercised the
default CV grid with three CV/final restarts and two native threads.
The TimeTree tests cover species-taxid resolution, missing child lineages,
ambiguous MRCA IDs, conflicting bounds, cached/offline retrieval and transport
failures. They also check that two-tip and other small clades remain eligible,
while size ranking and the query cap still apply. The workflow test switches from manual to recorded synthetic TimeTree
bounds without recomputing the species tree. The nwkit HTTP-client boundary
test requires the pinned `workflow/envs/timetree.yaml` installation; it replaces
only the HTTP response, so routine tests do not contact TimeTree. A separate
live API smoke test was run against four real species; its input topology and
branch lengths were supplied for integration testing, not inferred from BUSCO.
These tests do not establish accuracy or resource requirements for 5,586 species.
A separate tlight input check processed four real species and 20 shared BUSCOs:
all 80 CDS passed cdskit preparation. The previous QC retained 18 loci and
excluded two with eight informative sites each. Repeating trimAl with the
simplified QC retains all 20; the original 18 alignments/column maps are
unchanged. This check stopped after trimming; the original and updated results
are in `results/phylogeny_cdskit_trimal_smoke/summary.json` and
`results/phylogeny_qc_simplification_smoke/summary.json`, respectively.
