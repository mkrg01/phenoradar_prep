# Taxonomic review with MonoPhy

The `taxonomy_audit` target uses **MonoPhy 1.3.2** to compare an existing rooted
BUSCO **species tree** with registered NCBI taxonomy. It reports non-monophyletic
groups and MonoPhy's intruder/outlier tips, linked to all associated run IDs.
**Gene trees are not audit inputs and are not assessed.**

This replaces the former custom single-intruder detector entirely. Its
`min_reference_species` and `max_plot_species` settings are rejected with a
migration message. The old reference envelopes, per-gene reproduction scores
and candidate context plots are no longer produced. Regenerating an owned audit
directory replaces the complete report, including obsolete files.

Reports do not establish which sample is mislabeled. No metadata, inference
output or contrast pair is changed. After review, manual exclusion remains a
separate [`filter_species` export](species_filter.md), configured through the
top-level `exclude_species` list. Phenotyped trees and contrast results remain
available as provisional results.

## Configuration and execution

```yaml
phylogeny:
  species_sets: [all]
taxonomy_audit:
  enabled: false
  ranks: [family, subfamily, tribe, subtribe, genus]
  outlierlevel: 0.5
  collapse_monophyletic: true
  mem_gb: 8
```

The explicit target works with `enabled: false`:

```bash
./run_pipeline.sh --software-deployment-method conda \
  --configfile config/mydata.yaml config/phylogeny.local.yaml \
  --cores 1 --resources mem_gb=8 -- taxonomy_audit
```

The target follows `phylogeny.species_sets`: `[all]`, `[phenotyped]`, or both.
Reports go into the selected branch's `taxonomy_audit/` directory. The separate
representative tree under `contrast/phylogeny/` is not an audit target.
With `enabled: true`, reports are included in `phylogeny`, and in `all` when
`phylogeny.enabled: true`.

`ranks` may include other named NCBI ranks, such as order, class or subgenus.
Each requested rank is evaluated separately. MonoPhy itself can accept arbitrary
group memberships, but this workflow currently constructs columns from named
NCBI ranks only; it does not accept `no rank` as a single ambiguous grouping.

`outlierlevel` is MonoPhy's core-clade threshold, a number in `(0, 1]`. It is the
fraction of tips assigned to the focal group within a candidate clade, not a
confidence value, gene agreement fraction or fraction of that group's members
retained. The default is MonoPhy's 0.5. A higher threshold does not necessarily
mean fewer or more reliable flags. `collapse_monophyletic` affects display only;
there is no tip-count limit on assessment or plots.

The rule uses `workflow/envs/monophy.yaml`. Its adjacent post-deploy script
installs the unmodified official CRAN MonoPhy 1.3.2 source, verified by SHA-256.
For offline installation, set `MONOPHY_SOURCE_ARCHIVE` to that same archive.
No external reference sequence database is downloaded.

Completed species trees are reused under normal Snakemake dependency checks.
If they do not yet exist, the target schedules inference. Use `--dry-run` to
inspect upstream work before execution. The species-tree inference itself may
require gene trees; the audit does not read or re-analyze them.

For an archived result whose upstream workflow inputs are unavailable, activate
the MonoPhy environment and run the script directly:

```bash
python workflow/scripts/taxonomy_audit.py \
  --tree results/ANALYSIS/phylogeny/species_tree.nwk \
  --tree-qc results/ANALYSIS/phylogeny/species_tree.json \
  --samples results/ANALYSIS/metadata/samples.tsv \
  --taxonomy resources/taxonomy/taxa.sqlite \
  --outdir results/ANALYSIS/phylogeny/taxonomy_audit
```

For phenotyped trees, supply their selection manifest and species-tree/QC files.
The script requires exact agreement between tree tips and manifest species,
unique runs, consistent within-species metadata, and a rooted tree matching the
single outgroup recorded in its QC. It preserves the source root and does not
choose a new one. Input and code checksums and R/package versions are recorded.

## Interpretation

MonoPhy checks each registered group against the descendants of its most recent
common ancestor. If those descendants contain other groups, the group is
non-monophyletic. Its outlier procedure may select a smaller core clade using
`outlierlevel` and the numbers of focal-group tips in descendant branches.
Members outside the selected core are outliers; other groups' tips within it
are intruders. Ties and multifurcations can prevent the core search from resolving
which part to select. There is no branch-support threshold or calibrated
probability of mislabeling.

An event is relative to a **focal group**. A tip registered in group A can be an
outlier of A and an intruder of B simultaneously. Both events are retained in
`candidates.tsv`; MonoPhy's native `TipStates` instead gives intruder precedence.
`focal_taxon` is the group whose monophyly is being checked, not an inferred or
corrected identity for the tip.

Missing ranks are omitted separately **before** each MonoPhy assessment and
recorded as `missing_taxonomy_rank`. They are not treated as intruders or pooled
into an artificial unknown group. The pruned assessment tree is exported for
that rank. If fewer than two annotated tips remain, that rank is not passed to
MonoPhy. A one-tip group is `Monotypic`: it cannot be evaluated for its own
monophyly, although its tip can be reported as an intruder of another group.
Species-level assessment is consequently uninformative for within-species
monophyly when the input contains one tip per species.

Neither no flag nor monophyly verifies identity. Grouped mislabels, uneven
sampling, outdated taxonomy, genuinely non-monophyletic groups, hidden paralogy
and tree estimation errors can affect results. In particular, a majority of
incorrect labels can cause a correct tip to be flagged. The tool exposes
conflicts for review; it does not resolve that ambiguity. Several RNA-seq runs
may feed one species-level CDS set, so the report cannot identify which run is
responsible without separate sequence analysis.

## Outputs and figures

| File | Contents |
| --- | --- |
| `taxonomy.tsv` | Species-tree tip IDs and registered group columns; names include taxids to avoid ambiguous homonyms |
| `samples.tsv` | Original run manifest with review status and candidate event IDs |
| `taxon_results.tsv` | Each rank/group's monophyly, original MRCA node number, member and intruder/outlier counts |
| `rank_status.tsv` | Each species/rank's own-group status and separate intruder/outlier indicators |
| `candidates.tsv` | One native MonoPhy event per species/rank/focal group/role, with all associated runs |
| `events.tsv` | The same events before metadata enrichment |
| `group_members.tsv` | All focal-group members, its outliers, and other groups' intruder tips |
| `plot_members.tsv` | Every species' color, role symbol and display representative at each rank |
| `ranks/<rank>/monophy.rds` | Unmodified native MonoPhy result object, when assessment was possible |
| `ranks/<rank>/native_results.tsv`, `native_tip_states.tsv` | Native result and tip-state tables |
| `ranks/<rank>/assessment_tree.nwk` | Exact tree passed to MonoPhy, after omitting missing-rank tips; native MRCA numbers refer to this tree |
| `ranks/<rank>/tree.pdf`, `tree.svg`, `display_tree.nwk` | Publication figures in vector formats and corresponding display topology |
| `ranks/<rank>/figure.json` | Figure dimensions, font sizes and displayed group counts |
| `ranks/<rank>/not_assessed.txt` | Explanation when fewer than two annotated tips remain |
| `engine.json`, `R_session.txt`, `summary.json` | Engine settings, R/package versions, provenance and limitations |

`candidates.tsv` is an intruder/outlier list, not a complete list of all affected
species. Start with `taxon_results.tsv` to see every non-monophyletic group,
including those for which no outlier was selected. `samples.tsv` distinguishes
`review_flag`, `non_monophyletic_group`, `no_flag` and `not_assessable`.

Figures use native MonoPhy results and **ape-derived cladogram geometry**.
Each rank is exported as a vector PDF with embedded fonts and a vector SVG.
There are no titles, subtitles, footer notes or taxids in the figure itself.
Scientific names are set in italic at 9 pt; group annotations occupy a separate
aligned column at 8.1 pt. Page dimensions follow the measured label widths and
displayed tip count, keeping row spacing and font sizes consistent.

Registered groups share colors in both the tip markers and group annotations.
Small circles indicate unflagged tips; triangles indicate intruders, squares
outliers and diamonds both roles. An A tip within B keeps A's color. Missing-rank
tips have gray circles and blank group annotations. Role definitions belong in
the manuscript figure caption; `plot_members.tsv` retains every label, taxid,
color and symbol mapping. Internal branches remain neutral gray; no ancestral
color is interpreted as an inferred taxonomic assignment.

Each requested rank has separate PDF and SVG files. Unflagged monophyletic groups can be
folded to one labeled representative; non-monophyletic groups and all flagged
tips remain expanded. The fold is a display operation only and never reruns the
assessment. Species with missing rank annotations remain visible in gray on the
full-tree figure even though they were omitted from that rank's assessment.
`plot_members.tsv` records all folds. Figures are cladograms; branch lengths and
support labels should be inspected in the source or assessment Newick files.

## Validation and references

Tests invoke real MonoPhy and compare its saved object with a fresh direct
`AssessMonophyly` call. They cover role/run linkage, multiple grouped intruders,
missing ranks, singleton groups, parameter changes, obsolete-report replacement,
input validation and isolated Snakemake reuse with no gene-tree files present.
The independent manual species filter is also tested. These validate software
integration, not sensitivity or specificity for the study's samples.

- [MonoPhy paper](https://doi.org/10.7717/peerj-cs.56)
- [CRAN package](https://CRAN.R-project.org/package=MonoPhy)
- [Official source repository](https://github.com/oschwery/MonoPhy)
