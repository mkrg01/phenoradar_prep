# Inputs and configuration

[Back to README](../README.md)

Use this guide to adapt the workflow to a dataset. Run all commands from the repository root.

## Inputs

The default input layout is:

```text
metadata/
  metadata.tsv
transcriptome_assembly/
  multispecies_summary/
    busco_full_longest_cds.tsv
  longest_cds/
    {species}_longestCDS.fa.gz
  amalgkit_quant/
    {species}/{run}/{run}_abundance.tsv
```

| Input | Required columns or format |
| --- | --- |
| Sample metadata | Tab-separated; `scientific_name`, `run`, `taxid` |
| BUSCO summary | Tab-separated; `Species`, `busco_cds_single`, `busco_cds_duplicated`, `busco_cds_fragmented`, `busco_cds_missing`, `busco_cds_total` |
| Coding sequences | Gzip-compressed FASTA, one file per species |
| Abundance estimates | Tab-separated; `target_id`, `tpm`, one file per run |

Run IDs must be unique in the metadata, and species names must be unique in the
BUSCO summary. Multiple runs per species are supported. All runs for the same
species must agree on its taxid.

`{species}` is the scientific name with spaces replaced by underscores. Protein
filenames used by ODB-mapper additionally replace hyphens with underscores.
**Identifiers inside FASTA records are preserved.** Gene identifiers must be
unique across species and match the abundance table's `target_id` values. The
workflow rejects species-name collisions and checks that ODB queries belong to
the input species in their chunk.

Paths used by ODB-mapper must not contain spaces or shell metacharacters because
of limitations in the upstream shell scripts.

## Configuration

For a new dataset, copy the [default configuration](../config/config.yaml) and edit
it. If you already have a configuration, keep using it rather than overwriting it:

```bash
cp config/config.yaml config/mydata.yaml
```

Set the input paths, optional taxonomy source, analysis name, and appropriate OrthoDB
taxonomic node. **The default node `3193` is dataset-specific**; review it before
processing another dataset or downloading reference data.

`config/mydata.yaml`, `config/*.local.yaml`, and `config/pilot_species.txt` are
ignored by Git so that personal settings and species selections stay out of the
shared workflow. These files remain available locally.

| Setting | Purpose |
| --- | --- |
| `analysis` | Output directory name under `results/`, `work/`, and `logs/` |
| `inputs.*` | Metadata, BUSCO, CDS, and abundance paths |
| `translation.table` | Genetic code table used for CDS translation; default `1` |
| `taxonomy.source` | Optional existing SQLite database to copy when `resources/taxonomy/taxa.sqlite` is missing; default `null` downloads NCBI taxonomy |
| `selection.busco_threshold` | Minimum complete BUSCO fraction; default `0.5` |
| `selection.species_list` | Optional text file of species IDs, one per line |
| `selection.missing_taxonomy` | `error` or `allow` for unresolved taxids |
| `odb.existing_results` | Optional local import snapshot directory; `null` runs ODB normally |
| `odb.node` | Positive integer OrthoDB node; determines the reference path `resources/orthodb/v12_<node>` |
| `odb.chunk_size` | Maximum species per ODB chunk; default `50` |
| `odb.threads`, `odb.batch_size` | Worker and batch counts per chunk; defaults `16` and `64`; batch size must be at least the worker count |
| `odb.mem_gb` | Memory estimate per chunk in whole GB used to schedule ODB alongside other steps; default `192` |
| `odb.min_free_gb` | Minimum free disk space before mapping, in GiB; default `750` |
| `odb.reference_min_free_gb` | Minimum free disk space before reference preparation, in GiB; default `200` |
| `odb.allow_nonlocal` | Allow ODB work on filesystems outside the supported local types; default `false` |
| `odb.keep_work` | Whether to retain successful ODB work directories |
| `tpm.multimap` | Policy for genes assigned to multiple orthogroups; see [TPM interpretation](outputs.md#tpm-interpretation) |
| `alignment.enabled` | Include untrimmed, all-copy alignments of every mapped OG in the full workflow; default `false` |
| `alignment.threads`, `alignment.mem_gb` | Threads per OG alignment and decimal-GB memory budget per alignment-branch job; defaults `4` and `8` |
| `kegg.enabled` | Include KEGG outputs in the default full workflow; default `false` |
| `kegg.threads`, `kegg.mem_gb` | CPU and decimal-GB memory budgets per species; defaults `4` and `8` |
| `kegg.ambiguity` | `duplicate` (default) adds full TPM to each accepted KO; `drop` excludes multi-KO genes; `error` rejects quantified multi-KO genes. Annotations always retain candidates |

ODB chunks can run concurrently when their combined CPU and memory estimates fit
the workflow budget. The defaults use 16 workers and 192 GB for up to 50 species
per chunk. These are initial allowances to check against measured peak memory;
`config/pilot.yaml` retains smaller settings for its two-species chunks.

Configuration contains dataset paths, analysis choices, and resource budgets.
`inputs.species_trait` defaults to `species_trait/species_trait.tsv` and is the
sole source of phenotype annotations. It requires `species` and the column
selected by `phylogeny.trait` or `contrast.trait` (both default to `C4`). Spaces
in species names become underscores, matching the pipeline's species IDs; duplicate normalized names
are rejected. `C4` accepts `0`, `1`, or missing values. Metadata phenotype
columns are ignored. The file is required only by branches that use traits.
See [contrast pairs](contrast_pairs.md) for execution and outputs.

`phylogeny.species_sets` defaults to `[all]`. Use `[phenotyped]` for all
BUSCO-selected species with a nonmissing `phylogeny.trait`, or `[all, phenotyped]`
for both inference runs. They retain separate results in `phylogeny/` and
`phylogeny_phenotyped/`; changing only this list reuses completed outputs.
The phylogeny, preparation, calibration and dating targets all follow this list.
See [species sets](phylogeny.md#species-sets-and-reusable-outputs) for details.

Generated storage locations are fixed:

| Contents | Location |
| --- | --- |
| Analysis results, temporary work, logs | `results/<analysis>/`, `work/<analysis>/`, `logs/<analysis>/` |
| Taxonomy snapshot | `resources/taxonomy/taxa.sqlite` |
| OrthoDB reference | `resources/orthodb/v12_<node>/` |
| KOfam/KEGG reference | `resources/kegg/snapshot_v1/` |
| Completed KOfam/KEGG downloads | `resources/kegg/downloads/` |
| TimeTree response cache | `resources/timetree_cache/` |
| Launcher-managed Snakemake cache | `.cache/` |

Missing reference data are prepared automatically when a requested branch needs
them. Existing snapshots and cached responses are reused without automatic
updates. The KEGG reference is only required for KEGG targets or `kegg.enabled: true`.
See [reference management](references.md) for preparation and intentional refreshes.
The [OG alignment branch](alignments.md) reuses ODB mappings and translated
proteins; it requires no additional reference database.

Executable names and interpreters are fixed in the
workflow and supplied by each rule's Conda environment. ASTRAL-IV uses the
official build at `resources/phylogeny_tools/bin/astral4_int128`.

The workflow fixes OrthoDB to its supported version `v12`, standalone LSD2 to
one thread, and the TimeTree delay to one second between uncached requests.
These are not user configuration options. For phylogeny settings such as marker
selection, trimming, rooting and calibrations, see [the phylogeny guide](phylogeny.md).

Remove the following keys from older overrides: `tools`, `paths`,
`taxonomy.database`, `odb.reference_dir`, `kegg.reference_dir`, `odb.version`,
`kegg.command`, `phylogeny.famsa_command`, `phylogeny.trimal_command`,
`phylogeny.veryfasttree_command`, `phylogeny.astral_command`,
`phylogeny.dating.command`, `phylogeny.dating.threads`,
`phylogeny.dating.timetree.python`, and
`phylogeny.dating.timetree.request_delay_seconds`, and
`phylogeny.dating.timetree.cache_dir`. They are rejected with a
migration message. Execution records still retain the commands, tool hashes and
reference versions actually used; `run.json` records configuration and workflow
source checksums.

For the one-time migration of completed ODB results, `odb.existing_results` points
to a directory containing `annotations.tsv` and `snapshot.json` (schema version 1,
OrthoDB version/node, original protein checksums keyed by species, and the annotation
checksum). Preparing this snapshot is a local migration step, not a workflow target.
With this setting, `mapping` and the full workflow skip ODB reference preparation
and mapping. CDS translation still runs to verify exact agreement with the original
protein inputs. Missing species, changed proteins, incompatible version/node, or
modified annotations stop the import. Species subsets are supported; annotation
rows outside the selected FASTA IDs are excluded and counted in `merge_qc.json`.
The usual merged mappings and TPM outputs are produced. The explicit `references`
target still prepares an ODB reference if requested.

Memory settings use positive integers in decimal GB (1 GB = 1000 MB). For example,
`odb.mem_gb: 128` keeps the same memory estimate as the former `odb.mem_mb: 128000`.
Update older configuration files by replacing `odb.mem_mb` with `odb.mem_gb` and
dividing the value by 1000. See [memory units](running.md#memory-units) for Slurm and
direct execution.

Paths are relative to the repository root unless absolute. When run directly,
`run_pipeline.sh` changes to its own directory. With `sbatch`, submit from the
repository root or use `sbatch --chdir=/path/to/phenoradar_prep`; the script uses
the job's working directory because Slurm executes a temporary copy of it.
The launcher keeps Snakemake's cache in `.cache/` and sets `XDG_CACHE_HOME`
accordingly; `PHENORADAR_CACHE_DIR` is no longer used. The Snakemake executable
can be supplied with `SNAKEMAKE_BIN`.
Both execution modes print shell commands and rerun incomplete outputs automatically.

The workflow loads `config/config.yaml` and overlays files supplied with
`--configfile` in order. Command-line overrides such as
`--config analysis=my_experiment` take precedence. Use separate analysis names
when you want to keep results for comparison.

The default analysis name is `full`. [config/pilot.yaml](../config/pilot.yaml)
sets the name to `pilot`, restricts species using a user-provided list, and reduces
chunk sizes and resource estimates. The names determine output directories;
the configuration determines which data are processed. See the [pilot run guide](running.md#pilot-run).

When specifying targets such as `prepare` or `references`, place them after `--`
at the end of the command. `--configfile` accepts multiple filenames, so a target
immediately following a configuration filename would otherwise be read as another
configuration file. Keep all options before `--`.

### Species selection

Complete BUSCO fraction is `(single + duplicated) / total`. Species meeting or
exceeding the threshold are retained. An optional species list further restricts
that set; every requested species must exist and pass the threshold.

Species absent from the BUSCO table are excluded automatically. They remain in
`metadata_all.tsv` with `selected=False` and blank BUSCO values; `selection.json`
records their scientific names in `missing_busco_species` and their run count in
`missing_busco_runs`. They require no CDS, abundance, or taxonomy lookup and are
omitted from the BUSCO histogram. Explicit species lists still require every
requested species to have BUSCO data and pass the threshold.

Invalid counts in existing BUSCO rows, duplicate runs, or missing input files for
selected samples cause an error. Unresolved taxids for species with BUSCO data
cause an error by default; missing individual taxonomic ranks are allowed. Metadata columns such as
`exclusion` do not apply additional filters.

`metadata_all.tsv` retains all metadata rows with a `selected` flag. `samples.tsv`
records the selected species, run IDs, and exact CDS and abundance paths used by
downstream steps. BUSCO selection is by species; the histogram counts runs.
