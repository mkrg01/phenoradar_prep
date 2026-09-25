# Taxonomy check with MonoPhy

[Documentation](index.md) · [Species exclusion](species_filter.md)

`taxonomy_check` compares rooted BUSCO trees with NCBI taxonomy using MonoPhy.
It reports non-monophyletic groups and intruder/outlier tips linked to sample
runs. Flags support review; they never exclude species automatically.

## Configuration and execution

Set these options before [preparing analysis](datasets.md#run-an-analysis):

```yaml
phylogeny:
  trees: [all]
  taxonomy_check:
    enabled: true
    ranks: [family, subfamily, tribe, subtribe, genus]
    outlierlevel: 0.5
    collapse_monophyletic: true
```

```bash
./run_analysis.sh submit --analysis analyses/analysis001 --target taxonomy_check
```

The enabled flag is required and includes reports in `all`. The target follows
`phylogeny.trees`, including missing inference. Representative trees are unsupported.
Review jobs default to 1 CPU/8 GB; inference steps use their own resources.

Choose unique NCBI `ranks`. `outlierlevel` is the required focal-group fraction
inside a candidate core clade, in `(0, 1]`; it is not confidence.
`collapse_monophyletic` affects figures only.

## Interpretation

An **outlier** belongs to the focal group but falls outside its selected core
clade; an **intruder** belongs to another group but falls inside. A tip can have
both roles for different groups. `focal_taxon` names the assessed group, not a
corrected identity.

Missing ranks are omitted per assessment. Single-tip groups cannot be assessed
for their own monophyly, and some non-monophyletic groups yield no outliers.
There is no branch-support filter or probability of mislabeling.

Flags or their absence do not verify identity: sampling, taxonomy, paralogy,
and tree error can explain conflicts. Species-level results cannot identify a
responsible RNA-seq run without further sequence analysis.

## Outputs and figures

Under `results/<analysis>/phylogeny/<set>/taxonomy_check/`:

| File | Use |
| --- | --- |
| `taxon_results.tsv` | Start here: monophyly and intruder/outlier counts by rank/group |
| `candidates.tsv` | Tip events by focal group/role, linked to runs |
| `samples.tsv`, `rank_status.tsv` | Run-level review and species/rank roles |
| `taxonomy.tsv`, `group_members.tsv` | Registered taxonomy and assessed membership |
| `ranks/<rank>/` | Assessment tree, native MonoPhy result, PDF/SVG figures, and displayed membership |
| `summary.json`, `engine.json`, `R_session.txt` | Settings and provenance |

Figure colors denote registered groups: triangles mark intruders, squares
outliers, diamonds both, and gray tips missing ranks. Unflagged monophyletic
groups may be collapsed; use source/assessment trees for lengths and supports.

## References

[MonoPhy paper](https://doi.org/10.7717/peerj-cs.56) ·
[CRAN package](https://CRAN.R-project.org/package=MonoPhy) ·
[Source](https://github.com/oschwery/MonoPhy)
