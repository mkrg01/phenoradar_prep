# Development and tests

[Documentation](index.md) · [Releases](releases.md)

The entry point is [workflow/Snakefile](../workflow/Snakefile). Rules, processing
code, and environments live in `workflow/rules/`, `workflow/scripts/`, and
`workflow/envs/`; [layout.py](../workflow/scripts/layout.py) defines shared paths.

## Running tests

Run from the repository root with `pytest`, `pandas`, `ete4`, `matplotlib`,
`PyYAML`, `numpy`, `biopython`, and the pinned cdskit from
[phylogeny.yaml](../workflow/envs/phylogeny.yaml):

```bash
python -m pytest -q tests
```

[environment.yaml](../environment.yaml) provides Snakemake 9.8.0 and pytest for
CI/container generation; it does not include all test dependencies. Tool-dependent
tests skip when executables are unavailable. Override their paths as needed:

```bash
SNAKEMAKE_BIN=/path/to/snakemake SEQKIT_BIN=/path/to/seqkit \
  python -m pytest -q tests
```

| Tests | Additional tools |
| --- | --- |
| Core and KEGG integration | Snakemake, seqkit; ODB/KofamScan are test substitutes |
| Alignments, phylogeny, dating | FAMSA, trimAl, VeryFastTree, ASTRAL, LSD2; override with `FAMSA_BIN`, `TRIMAL_BIN`, `VERYFASTTREE_BIN`, `ASTRAL_BIN`, `LSD2_BIN` |
| Taxonomy audit | R/MonoPhy from `workflow/envs/monophy.yaml` |
| TimeTree client | Pinned nwkit from `workflow/envs/timetree.yaml` |

`PHYLOGENY_CONDA_PREFIX` enables tests with actual stage environments.
[Container CI](containers.md#build-and-publish-with-github-actions) runs
[smoke checks](../tests/container_smoke.py) with real tools;
[tests/container/Snakefile](../tests/container/Snakefile) checks Apptainer activation.

Release tests use temporary Git repositories and simulated GitHub/registry
responses. Slurm tests simulate allocations without submitting jobs. Use a real
dataset pilot to assess biological outputs and resource requirements.
