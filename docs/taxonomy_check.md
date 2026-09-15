# Taxonomy check with MonoPhy

[Documentation](index.md) · [Species exclusion](species_filter.md)

`taxonomy_check` compares rooted BUSCO trees with NCBI taxonomy using MonoPhy.
It reports non-monophyletic groups and intruder/outlier tips linked to sample
runs. Flags support review; they never exclude species automatically.

## Configuration and execution

```yaml
phylogeny:
  species_sets: [all]
taxonomy_check:
  enabled: false
  ranks: [family, subfamily, tribe, subtribe, genus]
  outlierlevel: 0.5
  collapse_monophyletic: true
```

```bash
./run_pipeline.sh --cores 1 --resources mem_gb=8 -- taxonomy_check
```

The target follows `phylogeny.species_sets` and can schedule missing inference:
check `--dry-run` and increase resources as needed. Set
`taxonomy_check.enabled: true` to attach reports to `phylogeny`.
Representative trees are not checked.

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

Under `results/<run_name>/phylogeny/<set>/taxonomy_check/`:

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

For archived results, [taxonomy_check.py](../workflow/scripts/taxonomy_check.py)
accepts source paths directly; see `--help`.

## References

[MonoPhy paper](https://doi.org/10.7717/peerj-cs.56) ·
[CRAN package](https://CRAN.R-project.org/package=MonoPhy) ·
[Source](https://github.com/oschwery/MonoPhy)
