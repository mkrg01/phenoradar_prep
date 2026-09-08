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

Set the input paths, taxonomy database, analysis name, and appropriate OrthoDB
taxonomic node. **The default node `3193` is dataset-specific**; review it before
processing another dataset or downloading reference data.

`config/mydata.yaml`, `config/*.local.yaml`, and `config/pilot_species.txt` are
ignored by Git so that personal settings and species selections stay out of the
shared workflow. These files remain available locally.

| Setting | Purpose |
| --- | --- |
| `analysis` | Output directory name under `results/`, `work/`, and `logs/` |
| `inputs.*` | Metadata, BUSCO, CDS, and abundance paths |
| `paths.*` | Root directories for results, temporary work, and logs |
| `tools.*` | Executable names or paths, plus an optional ODB environment prefix |
| `translation.table` | Genetic code table used for CDS translation; default `1` |
| `taxonomy.database` | Frozen ETE4 taxonomy SQLite database; see [reference setup](references.md#taxonomy-reference) |
| `selection.busco_threshold` | Minimum complete BUSCO fraction; default `0.5` |
| `selection.species_list` | Optional text file of species IDs, one per line |
| `selection.missing_taxonomy` | `error` or `allow` for unresolved taxids |
| `odb.version` | OrthoDB version; this workflow supports `v12` |
| `odb.node`, `odb.reference_dir` | OrthoDB node and reference snapshot directory |
| `odb.chunk_size` | Maximum species per ODB chunk; default `250` |
| `odb.threads`, `odb.batch_size` | Worker and batch counts; batch size must be at least the worker count |
| `odb.mem_gb` | Memory estimate in whole GB used to schedule ODB alongside other steps; default `128` |
| `odb.min_free_gb` | Minimum free disk space before mapping, in GiB; default `750` |
| `odb.reference_min_free_gb` | Minimum free disk space before reference preparation, in GiB; default `200` |
| `odb.allow_nonlocal` | Allow ODB work on filesystems outside the supported local types; default `false` |
| `odb.keep_work` | Whether to retain successful ODB work directories |
| `tpm.multimap` | Policy for genes assigned to multiple orthogroups; see [TPM interpretation](outputs.md#tpm-interpretation) |

Memory settings use positive integers in decimal GB (1 GB = 1000 MB). For example,
`odb.mem_gb: 128` keeps the same memory estimate as the former `odb.mem_mb: 128000`.
Update older configuration files by replacing `odb.mem_mb` with `odb.mem_gb` and
dividing the value by 1000. See [memory units](running.md#memory-units) for Slurm and
direct execution.

Paths are relative to the repository root unless absolute. When run directly,
`run_pipeline.sh` changes to its own directory. With `sbatch`, submit from the
repository root or use `sbatch --chdir=/path/to/phenoradar_prep`; the script uses
the job's working directory because Slurm executes a temporary copy of it.
The launcher keeps Snakemake's cache in `.cache/`. Override the cache location
with `PHENORADAR_CACHE_DIR`, or the Snakemake executable with `SNAKEMAKE_BIN`.
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

Missing BUSCO entries, invalid counts, duplicate runs, or missing input files for
selected samples cause an error. Unresolved taxids cause an error by default;
missing individual taxonomic ranks are allowed. Metadata columns such as
`exclusion` do not apply additional filters.

`metadata_all.tsv` retains all metadata rows with a `selected` flag. `samples.tsv`
records the selected species, run IDs, and exact CDS and abundance paths used by
downstream steps. BUSCO selection is by species; the histogram counts runs.
