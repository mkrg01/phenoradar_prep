"""Canonical result paths shared by producers and completed-result consumers."""

ORTHOGROUP_MAPPING = "orthogroups/mapping"
ORTHOGROUP_EXPRESSION = "orthogroups/expression"
ORTHOGROUP_ALIGNMENTS = "orthogroups/alignments"
PHYLOGENY_BRANCHES = {"all": "phylogeny/all", "phenotyped": "phylogeny/phenotyped"}
REPRESENTATIVES = "phylogeny/representatives"
ROOTING = f"{PHYLOGENY_BRANCHES['all']}/rooting"
CONTRAST_BRANCHES = frozenset(f"{branch}/contrast" for branch in
                             [*PHYLOGENY_BRANCHES.values(), REPRESENTATIVES])
