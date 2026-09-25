# Development and tests

[Documentation](index.md) · [Releases](releases.md)

The entry point is [workflow/Snakefile](../workflow/Snakefile). Rules, processing
code, and environments live in `workflow/rules/`, `workflow/scripts/`, and
`workflow/envs/`; [layout.py](../workflow/scripts/layout.py) defines shared paths.

## Running tests

Create the [test environment](../tests/environment.yaml) once. Run from the
repository root:

```bash
conda env create -n phenoradar-workflow -f tests/environment.yaml
conda activate phenoradar-workflow
python -m pytest -q -ra tests
```

If `phenoradar-workflow` already exists, update it first:

```bash
conda env update -n phenoradar-workflow -f tests/environment.yaml
```

[environment.yaml](../environment.yaml) is a minimal environment for Snakemake
and container generation. The full suite and CI use `tests/environment.yaml`,
which also provides ETE, plotting libraries, cdskit, and nwkit.
Tool-dependent tests skip when executables are unavailable; `-ra` lists the
reasons. Override their paths as needed:

```bash
SNAKEMAKE_BIN=/path/to/snakemake SEQKIT_BIN=/path/to/seqkit \
  python -m pytest -q tests
```

| Tests | Additional tools |
| --- | --- |
| Core and KEGG integration | Snakemake, seqkit; ODB/KofamScan are test substitutes |
| Alignments, phylogeny, dating | FAMSA, trimAl, VeryFastTree, ASTRAL, treePL; override with `FAMSA_BIN`, `TRIMAL_BIN`, `VERYFASTTREE_BIN`, `ASTRAL_BIN`, `TREEPL_BIN` |
| Taxonomy check | R/MonoPhy from `workflow/envs/monophy.yaml` |
| TimeTree client | Pinned nwkit from `workflow/envs/timetree.yaml` |

`PHYLOGENY_CONDA_PREFIX` enables tests with actual stage environments.
[GitHub Actions](../.github/workflows/container.yml) runs the same Python suite
and generates the Dockerfile on every branch push and pull request. The generated
file is saved as `container-recipe`, not tracked in Git. Manual runs and releases
build that recipe and run
[smoke checks](../tests/container_smoke.py) with real tools;
[tests/container/Snakefile](../tests/container/Snakefile) checks Apptainer activation.
For local image builds, first run `python workflow/scripts/generate_container.py`
in the pinned workflow environment; the generated Dockerfile is not committed.

Release tests use temporary Git repositories and simulated GitHub/registry
responses. Slurm tests simulate allocations without submitting jobs. Use a real
dataset pilot to assess biological outputs and resource requirements.
