# Taxonomic review with MonoPhy

[Documentation](index.md) · [Species exclusion](species_filter.md)

`taxonomy_audit` compares rooted BUSCO species trees with NCBI taxonomy using
MonoPhy. It reports non-monophyletic groups and intruder/outlier tips, linked to
sample run IDs. Flags support review and do not automatically exclude species.

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

```bash
./run_pipeline.sh --configfile config/mydata.yaml \
  --cores 1 --resources mem_gb=8 -- taxonomy_audit
```

The target follows `phylogeny.species_sets` and can schedule missing tree inference;
check `--dry-run` and increase resources if needed. Set `enabled: true` to attach
reports to `phylogeny`. Representative trees are not audited.

Choose named NCBI `ranks`. `outlierlevel` is the required focal-group fraction
within a candidate core clade, in `(0, 1]`; it is not confidence.
`collapse_monophyletic` changes figures only. The container supplies MonoPhy;
for offline native installation, `MONOPHY_SOURCE_ARCHIVE` can supply its pinned
archive from [monophy.yaml](../workflow/envs/monophy.yaml)'s installer.

For archived results, [taxonomy_audit.py](../workflow/scripts/taxonomy_audit.py)
also accepts tree, QC, sample-manifest, and taxonomy paths directly (`--help`).
Tree tips and the manifest must agree, including the recorded root/outgroup.

## Interpretation

An **outlier** is a focal group's member outside its selected core clade; an
**intruder** belongs to another group but lies inside that core. A tip can have
both roles for different focal groups. `focal_taxon` is the assessed group,
not a corrected identity. Some non-monophyletic groups yield no selected outlier.

Missing ranks are omitted separately for each assessment. A one-tip group cannot
be assessed for its own monophyly, but its tip can intrude into another group.
There is no branch-support filter or calibrated probability of mislabeling.

Neither monophyly nor absence of flags verifies identity. Sampling, taxonomy,
paralogy, and tree errors can produce conflicts. Species-level results cannot
identify which RNA-seq run is responsible without further sequence analysis.

## Outputs and figures

Reports are in `results/<run_name>/phylogeny/<set>/taxonomy_audit/`.

| File | Contents |
| --- | --- |
| `taxon_results.tsv` | Monophyly and intruder/outlier counts for every rank/group; start here |
| `candidates.tsv` | Tip events by focal group and role, linked to runs |
| `samples.tsv`, `rank_status.tsv` | Run-level review status and species/rank roles |
| `taxonomy.tsv`, `group_members.tsv` | Registered taxonomy and assessed group membership |
| `ranks/<rank>/assessment_tree.nwk`, `monophy.rds` | Exact assessed tree and native MonoPhy result |
| `ranks/<rank>/tree.pdf`, `tree.svg` | Figures; `plot_members.tsv` records displayed membership |
| `summary.json`, `engine.json`, `R_session.txt` | Settings and provenance |

Colors denote registered groups; triangles mark intruders, squares outliers, and
diamonds both. Unflagged monophyletic groups may be folded for display. Missing-rank
tips remain gray. Use source/assessment trees for branch lengths and support.

## References

- [MonoPhy paper](https://doi.org/10.7717/peerj-cs.56)
- [CRAN package](https://CRAN.R-project.org/package=MonoPhy)
- [Source](https://github.com/oschwery/MonoPhy)
