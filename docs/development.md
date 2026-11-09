# Development and tests

[Documentation](index.md) · [Releases](releases.md)

[workflow/Snakefile](../workflow/Snakefile) is the entry point. Rules, processing
code, and environments live in `workflow/rules/`, `workflow/scripts/`, and
`workflow/envs/`. [layout.py](../workflow/scripts/layout.py) defines shared paths.

## Running tests

From the repository root:

```bash
conda env create -n phenoradar-workflow -f tests/environment.yaml
conda activate phenoradar-workflow
python -m pytest -q -ra tests
```

For an existing environment, use
`conda env update -n phenoradar-workflow -f tests/environment.yaml`.
The test environment includes dependencies absent from the minimal host environment.
Tool-dependent tests skip when executables are unavailable; `-ra` lists the reasons.

| Tests | Additional tools or overrides |
| --- | --- |
| Core and KEGG | Snakemake, seqkit (`SNAKEMAKE_BIN`, `SEQKIT_BIN`); ODB/KofamScan use test substitutes |
| Alignments, phylogeny, dating | `FAMSA_BIN`, `TRIMAL_BIN`, `VERYFASTTREE_BIN`, `ASTRAL_BIN`, `LSD2_BIN` |
| Taxonomy check | R/MonoPhy from `workflow/envs/monophy.yaml` |
| TimeTree client | Pinned nwkit from `workflow/envs/timetree.yaml` |

`PHYLOGENY_CONDA_PREFIX` enables tests with actual stage environments. Scheduler
tests simulate Slurm/AGE services; AGE tests also run a real Snakemake DAG through
the shared launcher. Use a dataset pilot to assess biological outputs and resources.
For group/container changes, run `tests/slurm_busco_smoke.py` from an allocated
host environment with a cached SIF, marker FASTA and species manifest (see `--help`).
It checks identical alignments/column maps, tree inputs/tools/arguments, and completed-stage reuse.
Tree byte equality is reported separately; `--require-identical-trees` makes it mandatory.

## Container checks

[CI](../.github/workflows/container.yml) runs the Python suite and generates a
Dockerfile on branch pushes and pull requests. Manual runs and releases also
build the image and run [real-tool smoke checks](../tests/container_smoke.py).
Docker smoke checks use the caller's UID/GID and a writable temporary home;
Orthologer rejects root execution.
The [container test workflow](../tests/container/Snakefile) checks Apptainer activation.
ODB checks execute the native FASTA checker and import a small protein FASTA;
the same check must succeed before reference preparation or mapping chunks run.
The probe uses a temporary local API-version fixture. Reference downloads and
real mapping retain their normal API and pinned reference.
The ODB environment deployment script binds Orthologer Python entry points to
that environment's Python and records the repaired file hashes.

For a local image build, run `python workflow/scripts/generate_container.py` in
the pinned workflow environment. The generated Dockerfile is a build artifact
and is not committed.
