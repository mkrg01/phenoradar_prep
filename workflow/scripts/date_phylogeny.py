#!/usr/bin/env python3
"""Date the rooted CASTLES-II angiosperm tree with treePL."""
import argparse
import json
import math
import os
from pathlib import Path
import re
import subprocess
import tempfile
import time

from common import atomic_writer, file_record, now, read_tsv, write_json, write_tsv
from infer_phylogeny import executable, read_tree
from phylogeny_outgroup import outgroup_ids, validate_root

TREEPL_COMMIT = "f41af04ae7cc830deadbe83a1217ed9feca60c86"
TREEPL_ARCHIVE_SHA256 = "9d77514f03fa1583c1b207183e60fe11e2812218cbf83a77e37c643424d423d4"
ROOT_MAX_MA = 154
SMOOTHING = 10
SEED = 12345


def validate_dated(path, original, constraints):
    """Validate native six-decimal lengths without altering any branch length."""
    dated = read_tree(path, original.leaf_names())
    dated.dist = 0.0  # The root has no parent edge.
    intern = {}
    def codes(tree):
        result = {}
        for node in tree.traverse("postorder"):
            key = ("tip", node.name) if node.is_leaf else ("node", tuple(sorted(result[c] for c in node.children)))
            result[node] = intern.setdefault(key, len(intern))
        return result
    a, b = codes(original), codes(dated)
    if a[original] != b[dated] or len(a) != len(b):
        raise ValueError("treePL changed the rooted topology")
    names = {code: node.name for node, code in a.items()}
    for node, code in b.items():
        node.name = names[code]
    bounds = {row["node"]: row for row in constraints}
    lower, upper, counts, ages = {}, {}, {}, {}
    for node in dated.traverse("postorder"):
        if node.is_leaf:
            lo = hi = age = 0.0
            counts[node] = 1
        else:
            # Upstream prints each edge to six decimal places. Check whether
            # those intervals admit an ultrametric tree and the supplied bounds;
            # this is validation only, with no branch-length reconciliation.
            tolerance = {c: 0.500001e-6 + 2 * math.ulp(c.dist) for c in node.children}
            lo = max(lower[c] + max(0.0, c.dist - tolerance[c]) for c in node.children)
            hi = min(upper[c] + c.dist + tolerance[c] for c in node.children)
            counts[node] = sum(counts[c] for c in node.children)
            age = sum((ages[c] + c.dist) * counts[c] for c in node.children) / counts[node]
        if node.name in bounds:
            lo = max(lo, bounds[node.name].get("min_age_ma", 0.0))
            hi = min(hi, bounds[node.name].get("max_age_ma", math.inf))
        if lo > hi + 8 * math.ulp(max(1.0, abs(lo), abs(hi))):
            raise ValueError(f"treePL output violates calibration or ultrametricity beyond native rounding at {node.name}")
        lower[node], upper[node], ages[node] = lo, hi, age
    height = ages[dated]
    if not math.isfinite(height) or height <= 0:
        raise ValueError("treePL returned an invalid time-tree height")
    return dated, {n.name: age for n, age in ages.items() if not n.is_leaf}, height


def newick_text(tree):
    from ete4.parser.newick import make_parser
    return tree.write(parser=make_parser(1, dist="%.17g"), format_root_node=True) + "\n"


def treepl_build(command):
    """Accept only the checked, unmodified build installed by dating.post-deploy.sh."""
    binary = file_record(command)
    manifest = Path(command).resolve().parent.parent / "share/treepl/build.json"
    if not manifest.is_file():
        raise ValueError("treePL build provenance missing; run workflow/envs/dating.post-deploy.sh")
    build = json.loads(manifest.read_text())
    if (build.get("source_patch_applied") is not False or build.get("commit") != TREEPL_COMMIT
            or build.get("archive_sha256") != TREEPL_ARCHIVE_SHA256
            or build.get("executable_sha256") != binary["sha256"]):
        raise ValueError("treePL executable is not the verified unmodified build; recreate the dating environment")
    return build, binary


def tree_constraints(tree, path):
    rows = read_tsv(path)
    if not rows:
        raise ValueError("no applicable AngioCal fossil minimum; inspect angiocal/fossils.tsv")
    constraints = {}
    for row in rows:
        tips = [row["left_species"], row["right_species"]]
        if tips[0] == tips[1] or any(tip not in tree for tip in tips):
            raise ValueError(f"invalid AngioCal anchors: {tips}")
        node = tree.common_ancestor(tips)
        age = float(row["min_age_ma"])
        if node in constraints or not math.isfinite(age) or not 0 < age < ROOT_MAX_MA:
            raise ValueError(f"duplicate or incompatible AngioCal minimum: {tips}, {age} Ma")
        constraints[node] = dict(row, node=node.name, min_age_ma=age)

    # Native treePL rejects a stated ancestral minimum below a descendant minimum.
    # Omit only those redundant ancestral minima from the optimizer. Still check
    # every original fossil bound in the returned tree and retain all provenance.
    inherited, native_min = {}, {}
    for node in tree.traverse("postorder"):
        child_min = max((inherited[child] for child in node.children), default=0.0)
        minimum = constraints.get(node, {}).get("min_age_ma", 0.0)
        inherited[node] = max(child_min, minimum)
        if minimum and minimum >= child_min:
            native_min[node] = minimum

    if tree not in constraints:
        sides = sorted(min(child.leaf_names()) for child in tree.children)
        constraints[tree] = dict(node=tree.name, left_species=sides[0], right_species=sides[1])
    constraints[tree]["max_age_ma"] = ROOT_MAX_MA
    lines = []
    for node, row in constraints.items():
        if node in native_min or node.is_root:
            lines.append(f"mrca = {node.name} {row['left_species']} {row['right_species']}")
        if node in native_min:
            lines.append(f"min = {node.name} {native_min[node]:.17g}")
        if node.is_root:
            lines.append(f"max = {node.name} {ROOT_MAX_MA}")
    return list(constraints.values()), lines, len(rows) - len(native_min)


def date(tree, provenance, calibrations, outdir, command="treePL"):
    source, tree = tree, read_tree(tree)
    qc = json.loads(Path(provenance).read_text())
    if qc.get("branch_length_unit") != "substitutions_per_site" or not qc.get("outgroup"):
        raise ValueError("dating requires a rooted tree with substitution-per-site lengths")
    outgroup = outgroup_ids(qc["outgroup"])
    validate_root(tree, outgroup)
    if any(len(node.children) != 2 for node in tree.traverse() if not node.is_leaf):
        raise ValueError("treePL dating requires a rooted bifurcating tree")
    tips = set(tree.leaf_names())
    if any(not re.fullmatch(r"[A-Za-z0-9_.-]+", tip) for tip in tips):
        raise ValueError("treePL tip labels must contain only letters, digits, underscores, dots or hyphens")
    for i, node in enumerate(tree.traverse(), 1):
        if not node.is_leaf:
            node.name = f"DatingNode{i}"
            while node.name in tips:
                node.name = "_" + node.name
    sites = qc.get("total_gene_sites")
    if type(sites) is not int or not 0 < sites < 2**31:
        raise ValueError("treePL requires total_gene_sites as a positive 32-bit integer")
    if not any(node.dist > 0 for node in tree.traverse() if not node.is_root):
        raise ValueError("treePL requires positive substitution lengths")
    constraints, native_bounds, redundant = tree_constraints(tree, calibrations)
    command = executable(command)
    build, binary = treepl_build(command)
    config = "\n".join([
        "treefile = input.nwk", "outfile = dated.nwk", f"numsites = {sites}",
        "nthreads = 1", f"seed = {SEED}", f"smooth = {SMOOTHING}",
        "opt = 1", "optad = 1", *native_bounds,
    ]) + "\n"
    outdir = Path(outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    start = time.monotonic()
    with tempfile.TemporaryDirectory(prefix=".treepl-", dir=outdir) as tmp:
        work = Path(tmp)
        (work / "input.nwk").write_text(newick_text(tree))
        (work / "treepl.conf").write_text(config)
        argv = [command, "treepl.conf"]
        native_tree = work / "dated.nwk"
        with (work / "treepl.log").open("w") as log:
            result = subprocess.run(argv, cwd=work, stdout=log, stderr=subprocess.STDOUT,
                                    env={**os.environ, "OMP_NUM_THREADS": "1", "OMP_THREAD_LIMIT": "1"})
        report = (work / "treepl.log").read_text()
        # Some native error paths exit 0 without producing a completed fit.
        fits = re.findall(r"^after opt calc:\s+(\S+)", report, flags=re.MULTILINE)
        if result.returncode or not native_tree.is_file() or len(fits) != 1:
            raise ValueError(f"treePL did not complete dating (exit {result.returncode}):\n{report[-5000:]}")
        objective = float(fits[0])
        if not math.isfinite(objective) or not 0 <= objective < 1e15:
            raise ValueError(f"treePL returned an invalid objective: {objective}")
        dated, ages, height = validate_dated(native_tree, tree, constraints)
        validate_root(dated, outgroup)
        with atomic_writer(outdir / "treepl.log") as handle:
            handle.write(report)
    with atomic_writer(outdir / "treepl.conf") as handle:
        handle.write(config)
    with atomic_writer(outdir / "species_tree.dated.nwk") as handle:
        handle.write(newick_text(dated))
    write_tsv(outdir / "node_ages.tsv", ["node", "age_ma"],
              [{"node": node, "age_ma": age} for node, age in ages.items()])
    write_json(outdir / "provenance.json", {
        "created_at": now(), "tree": file_record(source), "source_provenance": file_record(provenance),
        "calibrations": file_record(calibrations), "build": build, "executable": binary,
        "method": "treePL penalized likelihood", "branch_length_unit": "million_years",
        "calibration_source": "AngioCal v1.0", "root_max_age_ma": ROOT_MAX_MA,
        "root_max_source": "Zuntini et al. (2024) young-tree crown-angiosperm upper bound; https://doi.org/10.1038/s41586-024-07324-0",
        "root_age_ma": height, "root_at_upper_bound": math.isclose(height, ROOT_MAX_MA, abs_tol=1e-4, rel_tol=0),
        "smoothing": SMOOTHING, "seed": SEED, "numsites": sites, "objective": objective,
        "tip_count": len(tips), "elapsed_seconds": time.monotonic() - start,
        "topology_preserved": True, "outgroup": outgroup, "redundant_fossil_minima": redundant,
        "zero_internal_edges": sum(node.dist == 0 for node in dated.traverse()
                                   if not node.is_root and not node.is_leaf),
        "optimizer_line_search_failures": report.count("Linear search failed"),
    })
    print(f"treePL: {len(tips)} tips, root {height:.6f} Ma", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["tree", "provenance", "calibrations", "outdir"]:
        parser.add_argument("--" + name, required=True)
    parser.add_argument("--command", default="treePL")
    date(**vars(parser.parse_args()))
