# Development and tests

[Documentation](index.md)

The workflow entry point is [workflow/Snakefile](../workflow/Snakefile).
Stage rules live in `workflow/rules/`, processing code in `workflow/scripts/`,
and software definitions in `workflow/envs/`. Shared producer/export paths are
defined in [layout.py](../workflow/scripts/layout.py).

## Running tests

Run tests from the repository root in a Python environment containing `pytest`,
`pandas`, `ete4`, `matplotlib`, `PyYAML`, `numpy`, `biopython`, and the pinned
cdskit from [phylogeny.yaml](../workflow/envs/phylogeny.yaml). The launcher
environment alone does not contain all test dependencies.

```bash
python -m pytest -q tests
```

Tool-dependent tests use installed executables or skip when unavailable. Supply
paths when needed:

```bash
SNAKEMAKE_BIN=/path/to/snakemake \
SEQKIT_BIN=/path/to/seqkit \
  /path/to/test-environment/bin/python -m pytest -q tests
```

| Tests | Additional software |
| --- | --- |
| Core workflow integration | Snakemake and seqkit; ODB is a test substitute |
| KEGG integration | Snakemake and seqkit; KofamScan is a test substitute |
| OG alignments | FAMSA (`FAMSA_BIN`); workflow tests also need Snakemake/seqkit |
| Phylogeny and dating | FAMSA, trimAl (`TRIMAL_BIN`), VeryFastTree (`VERYFASTTREE_BIN`), ASTRAL (`ASTRAL_BIN`), LSD2 (`LSD2_BIN`) |
| Taxonomy audit | `Rscript` with MonoPhy and dependencies from `workflow/envs/monophy.yaml` |
| TimeTree HTTP-client boundary | Pinned nwkit from `workflow/envs/timetree.yaml` |

The phylogeny integration uses the verified ASTRAL build at its fixed project
location. `PHYLOGENY_CONDA_PREFIX` enables testing with actual stage Conda
environments. `LSD2_SOURCE_ARCHIVE` can supply the verified LSD2 source offline.
The executable variables above select test tools; production commands are fixed.

The TimeTree adapter calls `nwkit.mcmctree._fetch_timetree_url` and validates
the pinned client-module hash. Changes to that private API require adapter
review; its boundary test substitutes only the HTTP response.

## Coverage

| Area | Main checks |
| --- | --- |
| Inputs and taxonomy | Species/run identity, BUSCO selection, missing data, local snapshots, download/build failure recovery |
| ODB and TPM | Output formats, protein/reference checksums, ambiguous mappings, normalization, concurrent chunks, interrupted jobs, abundance-only updates |
| KEGG | Reference inventories and download recovery, hit acceptance, ambiguity policies, gene-level coverage, overlapping KO sums, zero versus missing observations |
| OG alignments | All copies, repeated assignments, singleton OGs, species-encoded IDs, residue preservation, stale OG removal, per-OG recovery |
| BUSCO phylogeny | Coverage ranking, sequence QC, cdskit API/CLI agreement, trimAl columns, original X retention, missing-species rejection, independent species sets and roots |
| Dating | Calibration MRCA/bounds, all three variance settings, zero branches, finite native output, rooted topology, time units, rounding, ultrametricity, dating-only recovery |
| TimeTree | Taxid resolution, child-lineage coverage, repeated MRCA IDs, small clades, query caps, conflicting bounds, cached/offline responses, transport failures |
| Contrast pairs | Representative maps, deterministic ties, unknown traits, unresolved/empty pairs, figures, direct molecular membership, pair reformation after exclusions |
| Taxonomic review | Real MonoPhy agreement, intruder/outlier roles, run linkage, missing ranks, singleton groups, report replacement |
| Filtering and collection | Completed-result operation without raw inputs, identity/checksum validation, retained values and alignment columns, pruned path lengths, reversible exclusions, failed collection preserving prior output |
| Launcher | Direct/Slurm argument handling, spooled script paths, resource budgets, exit status, real local Snakemake execution |

Integration fixtures create isolated projects, expose tools under the workflow's
fixed command names, and exercise unchanged reruns and targeted recovery. Tests
using synthetic reference fixtures and substituted HTTP responses do not need
full database downloads or live TimeTree queries. Slurm tests simulate an
allocation; they do not submit cluster jobs.

## Validation limits

These checks establish processing behavior and integration. Test substitutes
cannot validate real ODB or KO assignments, synthetic phylogenies do not measure
biological accuracy, and small inputs do not establish full-scale resource use.
Run a real pilot for the intended dataset.

Historical real-input checks are recorded in [dataset checks](notes/validation.md).
The [method comparison](notes/phylogeny_comparison.md) and
[dating-tool evaluation](notes/dating_evaluation.md) retain methodological review
and dated measurements separately from the test instructions.
