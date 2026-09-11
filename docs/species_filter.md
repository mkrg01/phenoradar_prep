# Manual species exclusion from completed outputs

Use a **top-level YAML list** of exact `species` IDs from the original
`metadata/samples.tsv`:

```yaml
exclude_species:
  - Lespedeza_davurica
  - Cleistogenes_squarrosa
  - Phragmites_karka
```

These names are an example, not an enabled default. The default is
`exclude_species: []`. There is no `species_filter:` configuration wrapper,
no automatic use of taxonomy-audit candidates, and no run-specific exclusion
setting. All runs and gene copies belonging to a listed species are excluded.
Unknown IDs, duplicates, invalid identifiers, and removal of every selected
species fail. IDs use underscores in place of spaces; existing hyphens remain.

## Execute after the source analysis

Save the list in your existing local configuration or a separate local override,
for example `config/exclusions.local.yaml`, then run:

```bash
./run_pipeline.sh --software-deployment-method conda \
  --configfile config/mydata.yaml config/exclusions.local.yaml \
  --cores 1 --resources mem_gb=8 -- filter_species
```

`filter_species` is an explicit postprocessing target. It is not automatically
included in `all`, `phylogeny`, or `contrast_pairs`. It reads completed files
under `results/<analysis>/` and writes `results/<analysis>/filtered/`.
The original selected-sample manifest must exist. A dry run is available with
`--dry-run` before `-- filter_species`.

Only completed output groups are exported. A branch with no final outputs is
`absent`; a partially completed group is `incomplete`, with missing filenames
listed in the new manifest. Neither starts upstream jobs. Missing files claimed
by a completed alignment inventory, or inconsistent identities/checksums within
completed inputs, fail rather than produce a partial success.

The rule uses absolute paths to completed snapshot files as external inputs.
The original relative-path producer rules are not dependencies of this target:
raw RNA-seq, CDS, BUSCO tables, external reference preparation, mapping and
phylogeny tools do not need to be available. The analysis Conda environment
supplies the Python/ETE4 code used to prune trees.

Changing the list regenerates only the export. Changed source files or a newly
completed optional branch invalidate it on the next invocation. Unchanged
inputs/list reuse the existing export. Every export starts from the original
analysis, so removing an ID from the exclusion list restores that species.
Using a previous filtered export as the source is rejected.

For archived results outside the ordinary workflow layout:

```bash
python workflow/scripts/filter_species.py \
  --source /path/to/completed/analysis \
  --exclude-species '["Lespedeza_davurica"]' \
  --traits species_trait/species_trait.tsv \
  --outdir /path/to/new/filtered
```

`--traits` and `--outdir` are optional; the latter defaults to `<source>/filtered`.
The standalone script validates and stages a complete export before replacing
its own previous directory. Original inputs and unrelated output directories
are never overwritten. Normal Snakemake failed-output cleanup still applies
when using the workflow target.

## Exported dataset

The curated dataset is a downstream input bundle, not a second independent
execution of the original workflow. `manifest.json` is its provenance record;
old selection counts, inference QC and execution claims are not copied under
the guise of new calculations.

| Output under `filtered/` | Behavior |
| --- | --- |
| `metadata/samples.tsv`, `metadata/species.txt` | Remaining selected samples/species; original sample columns and run IDs retained |
| `metadata/metadata_all.tsv`, `metadata/metadata_high_busco.tsv` when available | Only remaining active sample rows; pre-BUSCO excluded rows are not reintroduced |
| `metadata/species_trait.tsv` when available | Remaining species' original phenotype columns, preserving zero and missing values; spaces in species IDs normalized |
| `excluded_samples.tsv` | Removed species/run identities and taxids for review |
| `proteins/` when complete | Symlinks to original files for remaining species, using the manifest's `odb_species` filename mapping |
| `odb/merged/mappings.sqlite`, `gene_orthogroups.tsv`, `filter_qc.json` | Subset of the original gene ownership/mapping database, including all retained copies and ambiguous gene/OG assignments |
| `tpm/*.tsv` | Original long/wide tables and run QC with excluded runs removed |
| `kegg/*.tsv` | Filtered gene/KO membership, expression, support and run QC; feature descriptions retained |
| `alignments/*.faa`, `members.tsv`, `filter_qc.json` | Retained gene rows with original full headers, residues and site coordinates; empty OGs omitted and listed |
| `phylogeny/species_tree.pruned.nwk` | All-species tree with excluded tips removed, when at least two tips remain |
| `phylogeny/gene_trees.pruned.nwk`, `gene_trees/*.pruned.nwk`, `gene_trees.tsv` | Per-marker derivatives and a marker/retention index; trees with fewer than two tips omitted |
| `phylogeny/alignments/`, including `raw/` when present | Available retained-marker alignments with excluded species rows removed; saved column maps unchanged |
| `phylogeny/species_coverage.tsv`, `pruning.json` | Recounted retained-tree coverage, root status and explicit pruning limitations |
| `phylogeny/dating/species_tree.dated.pruned.nwk` when available | Pruned source time tree, without refitting ages or calibration constraints |
| `manifest.json` | Exclusion list, retained identities, before/after counts, exported/skipped branches and input/output/code checksums |

Gene ownership comes from database or membership rows, never gene-name prefix
parsing. Where ODB mappings exist, alignment and KO gene ownership must agree
with that database. Original per-run outputs, chunk results, logs and external
CDS/abundance inputs remain source caches; they are not duplicated. Paths in
the selected-sample manifest still point to the original inputs.

The full source phenotype table is unchanged. `phylogeny_phenotyped/` and
`contrast/` remain provisional original results and are not copied into the
curated bundle. Taxonomy-audit reports stay with the source analysis as the
record of why a species was considered for exclusion. Future contrast-pair
selection should explicitly consume the curated full-species dataset.

## Numerical and phylogenetic meaning

Expression values are copied as strings without reaggregation or normalization.
OG TPM is normalized within each run, so removing other runs does not change
the retained runs' denominators. KO zeros and unavailable values remain distinct.
OG/KO feature columns retain the original axes, including all-zero or empty
columns; this stage does not perform feature selection. An OG can therefore
remain a table column after its last alignment sequence has been removed.
Such empty alignments are listed in `alignments/filter_qc.json`.

Alignment row removal retains every column, including columns that become all
gap. It is not a new alignment or a new trimAl pass. Marker selection is not
rerun. Original alignments remain available if a later analysis needs to revisit
alignment or site selection after exclusion.

Pruned trees retain path lengths by summing contracted edges. They are **not
new species-tree, gene-tree or dating estimates**: the original inference can
still reflect the excluded sequences. Internal labels/supports are omitted
because supports on the contracted tree have not been recomputed. Original
support values remain in the source trees. If the original outgroup is removed,
`pruning.json` marks the root as requiring review; no new biological root is
silently inferred. `gene_trees.tsv:at_least_four_tips` records tip-count eligibility
only; it does not verify suitability for an inference program's other settings.
Trees with two or three tips are retained as derivatives. Old node IDs, age tables, calibration tables
and inference-QC files are not relabeled as filtered calculations.

The saved BUSCO alignments permit later tree inference without redoing RNA-seq
assembly, translation or ODB mapping. This export does not automatically run
that inference or finalize new contrast pairs.

## Storage and verification

Large tabular inputs are streamed or joined through SQLite. Filtered tables and
the mapping database require new disk space; this is data export rather than
sequence reanalysis. Protein files use symlinks to avoid duplicating the large
FASTA collection. Keep original results available and treat linked protein
files as shared input data; editing through a link would edit the original.
Input/output checksum verification still reads those files.

Tests cover replicated runs, unrelated gene IDs, multiple copies and ambiguous
assignments, numeric-string/zero/missing-value preservation, alignment coordinates,
empty OGs, path-length preservation, missing outgroups, small trees, exclusion
reversal, incomplete branches and invalid/stale input rejection. A full-Snakefile
test runs the export without raw inputs or inference/reference tools and verifies
that changing the list/source regenerates only the export.
