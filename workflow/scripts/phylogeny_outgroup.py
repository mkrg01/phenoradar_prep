"""Outgroup sample sets shared by inference and downstream tree consumers."""
from pathlib import Path


def outgroup_ids(value):
    """Normalize current lists and historical single-sample QC records."""
    if isinstance(value, str):
        value = [value]
    if (not isinstance(value, (list, tuple)) or not value
            or any(not isinstance(n, str) or not n or any(c.isspace() for c in n) for n in value)
            or len(set(value)) != len(value)):
        raise ValueError("outgroup must contain distinct nonempty sample IDs")
    return sorted(value)


def read_outgroup(path):
    return outgroup_ids(Path(path).read_text().splitlines())


def validate_outgroup(names, outgroup):
    outgroup = set(outgroup_ids(outgroup))
    names = set(names)
    if not outgroup <= names:
        raise ValueError("outgroup is absent from the selected sample manifest: " + ", ".join(sorted(outgroup - names)))
    if outgroup == names:
        raise ValueError("outgroup leaves no ingroup samples")
    return outgroup


def validate_root(tree, outgroup):
    outgroup = validate_outgroup(tree.leaf_names(), outgroup)
    if len(tree.children) != 2 or not any(set(c.leaf_names()) == outgroup for c in tree.children):
        raise ValueError("species-tree root does not match the selected outgroup clade")


def root_on_outgroup(tree, outgroup):
    """Orient an existing split without pruning tips or changing unrooted edges.

    The midpoint used for Newick serialization is not an estimated time split;
    LSD2 re-estimates the position along this edge during dating.
    """
    wanted = validate_outgroup(tree.leaf_names(), outgroup)
    if len(tree.children) == 2 and any(set(c.leaf_names()) == wanted for c in tree.children):
        return tree
    sizes, hits, candidates = {}, {}, []
    total = len(list(tree.leaf_names()))
    for node in tree.traverse("postorder"):
        sizes[node] = 1 if node.is_leaf else sum(sizes[c] for c in node.children)
        hits[node] = int(node.name in wanted) if node.is_leaf else sum(hits[c] for c in node.children)
        if not node.is_root and ((sizes[node] == hits[node] == len(wanted))
                                or (hits[node] == 0 and sizes[node] == total - len(wanted))):
            candidates.append(node)
    if not candidates:
        raise ValueError("selected basal group is not monophyletic: no edge separates all outgroup samples")
    tree.set_outgroup(candidates[0])
    validate_root(tree, outgroup)
    return tree
