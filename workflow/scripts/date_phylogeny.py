#!/usr/bin/env python3
"""Date CASTLES substitution lengths with LSD2 on the selected basal edge."""
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
from phylogeny_outgroup import outgroup_ids, root_on_outgroup, validate_root

DEFAULT_SETTINGS = {"variance": 1, "variance_parameter": None}
LSD2_COMMIT = "c61110f3a4fa05325b45c97b2134792ff9d55d4c"
LSD2_ARCHIVE_SHA256 = "9bbeaa0f8f35783c1d8dec74df6c93a804dbca808fa04484f9123de4e7258b53"


def validate_settings(settings):
    if not isinstance(settings, dict):
        raise ValueError("LSD2 settings must be a mapping")
    unknown = set(settings) - set(DEFAULT_SETTINGS)
    if unknown:
        raise ValueError(f"unknown LSD2 settings: {sorted(unknown)}")
    cfg = {**DEFAULT_SETTINGS, **settings}
    if type(cfg["variance"]) is not int or cfg["variance"] not in (0, 1, 2):
        raise ValueError("LSD2 variance must be 0, 1 or 2")
    b = cfg["variance_parameter"]
    if b is not None and (type(b) not in (int, float) or not math.isfinite(b) or not 0 < b <= 1):
        raise ValueError("LSD2 variance_parameter must be null or finite in (0, 1]")
    if cfg["variance"] == 0 and b is not None:
        raise ValueError("LSD2 variance_parameter requires variance 1 or 2")
    return cfg


def calibration_rows(tree, path):
    names = set(tree.leaf_names())
    constraints, used = [], set()
    for row in read_tsv(path):
        if not {"taxa", "min_age_ma", "max_age_ma", "source"} <= row.keys():
            raise ValueError("calibrations require taxa, min_age_ma, max_age_ma, source columns")
        taxa = row["taxa"].split(",")
        if len(set(taxa)) != len(taxa) or len(taxa) < 2 or set(taxa) - names:
            raise ValueError(f"calibration taxa must name at least two distinct tree tips: {row['taxa']}")
        minimum, maximum = float(row["min_age_ma"]), float(row["max_age_ma"])
        if not (math.isfinite(minimum) and math.isfinite(maximum) and 0 < minimum <= maximum):
            raise ValueError("calibration ages must be finite positive Ma bounds with min <= max")
        if not row["source"].strip():
            raise ValueError("every calibration must record its source")
        node = tree.common_ancestor(taxa)
        if node.name in used:
            raise ValueError(f"multiple calibrations resolve to the same MRCA: {row['taxa']}")
        used.add(node.name)
        constraints.append({**row, "node": node.name, "min_age_ma": minimum, "max_age_ma": maximum})
    if not constraints:
        raise ValueError("absolute dating requires at least one explicit calibration; no arbitrary root age is supplied")
    bounds = {c["node"]: c for c in constraints}
    for node in tree.traverse("postorder"):
        lower = max((child.props["age_lower"] for child in node.children), default=0)
        bound = bounds.get(node.name)
        if bound:
            lower = max(lower, bound["min_age_ma"])
            if lower > bound["max_age_ma"]:
                raise ValueError(f"inconsistent ancestor/descendant calibration bounds at {node.name}")
        node.add_prop("age_lower", lower)
    return constraints


def native_rounding(value):
    """Half a unit of the last digit in native default C++ stream precision."""
    return (0.500001 * 10 ** (math.floor(math.log10(abs(value))) - 5)
            if value else 0.0) + 2 * math.ulp(value)


def validate_dated(path, original, constraints):
    """Validate native six-significant-digit lengths without altering any branch length."""
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
        raise ValueError("LSD2 changed the rooted topology")
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
            # Upstream prints each edge to six significant digits. Check whether
            # those intervals admit an ultrametric tree and the supplied bounds;
            # this is validation only, with no branch-length reconciliation.
            tolerance = {c: native_rounding(c.dist) for c in node.children}
            lo = max(lower[c] + max(0.0, c.dist - tolerance[c]) for c in node.children)
            hi = min(upper[c] + c.dist + tolerance[c] for c in node.children)
            counts[node] = sum(counts[c] for c in node.children)
            age = sum((ages[c] + c.dist) * counts[c] for c in node.children) / counts[node]
        if node.name in bounds:
            lo = max(lo, bounds[node.name]["min_age_ma"])
            hi = min(hi, bounds[node.name]["max_age_ma"])
        if lo > hi + 8 * math.ulp(max(1.0, abs(lo), abs(hi))):
            raise ValueError(f"LSD2 output violates calibration or ultrametricity beyond native rounding at {node.name}")
        lower[node], upper[node], ages[node] = lo, hi, age
    height = ages[dated]
    if not math.isfinite(height) or height <= 0:
        raise ValueError("LSD2 returned an invalid time-tree height")
    return dated, {n.name: age for n, age in ages.items() if not n.is_leaf}, height


def newick_text(tree):
    from ete4.parser.newick import make_parser
    return tree.write(parser=make_parser(1, dist="%.17g"), format_root_node=True) + "\n"


def lsd2_build(command):
    binary = file_record(command)
    manifest = Path(command).resolve().parent.parent / "share/lsd2/build.json"
    if not manifest.is_file():
        raise ValueError("LSD2 build provenance missing; run workflow/envs/dating.post-deploy.sh")
    build = json.loads(manifest.read_text())
    if (build.get("source_patch_applied") is not False or build.get("commit") != LSD2_COMMIT
            or build.get("archive_sha256") != LSD2_ARCHIVE_SHA256
            or build.get("executable_sha256") != binary["sha256"]):
        raise ValueError("LSD2 executable is not the verified unmodified build; recreate the dating environment")
    return build, binary


def dates_text(constraints):
    # LSD2 dates increase towards the present; positive ages in Ma run backwards.
    lines = []
    for row in constraints:
        young, old = row["min_age_ma"], row["max_age_ma"]
        bound = f"{-young:.17g}" if young == old else f"b({-old:.17g},{-young:.17g})"
        lines.append(f"mrca({row['taxa']}) {bound}")
    return str(len(lines)) + "\n" + "\n".join(lines) + "\n"


def native_result(path, variance):
    text = Path(path).read_text()
    matches = re.findall(r"^ rate (\S+), tMRCA (\S+) , objective function (\S+)$", text, re.M)
    if len(matches) != (2 if variance == 2 else 1) or "TOTAL ELAPSED TIME:" not in text:
        raise ValueError(f"LSD2 final optimization is incomplete; see {path}")
    rate, root, objective = matches[-1]
    rates, roots = [float(v) for v in rate.split(":")], [float(v) for v in root.split(":")]
    objective = float(objective)
    if (len(rates) not in (1, 2) or len(roots) != len(rates)
            or not all(math.isfinite(v) and v > 0 for v in rates)
            or not all(math.isfinite(v) and v < 0 for v in roots)
            or not math.isfinite(objective) or objective < 0):
        raise ValueError(f"LSD2 returned an invalid or unbounded rate/date solution; see {path}")
    warnings = re.findall(r"\*WARNINGS:\s*(.*?)\*RESULTS:", text, re.S)
    return {"objective": objective, "rate_substitutions_per_site_per_ma": rates,
            "root_age_solution_bounds_ma": sorted(-v for v in roots),
            "unique_time_scale": len(roots) == 1,
            "diagnostics": [w.strip() for w in warnings if w.strip()]}


def extract_time_tree(nexus):
    # The .nwk output is in substitutions/site. Only .date.nexus has time edges.
    records = re.findall(r"^tree\s+\S+\s*=\s*(.+;)\s*$", Path(nexus).read_text(), re.M | re.I)
    if len(records) != 1:
        raise ValueError("LSD2 must return exactly one native time tree")
    return re.sub(r"\[[^\]]*\]", "", records[0]) + "\n"


def normalize_numerical_zeros(text, root_age):
    """Clamp only tiny negative solver residuals at zero, keeping native output.

    Temporal constraints allow equality. Floating-point cancellation can emit
    negative lengths of a few ulps despite a zero lower bound. Larger negatives
    are errors, not branches to repair. Report every changed edge.
    """
    from ete4 import Tree
    tree = Tree(text.strip(), parser=1)
    tolerance = 64 * math.ulp(max(1.0, root_age))
    adjustments = []
    for number, node in enumerate(tree.traverse("preorder"), 1):
        if not node.is_root and node.dist is not None and node.dist < 0:
            if not math.isfinite(node.dist) or node.dist < -tolerance:
                raise ValueError("LSD2 returned a negative time branch beyond floating-point roundoff")
            adjustments.append({"preorder_index": number, "native_age_length_ma": node.dist,
                                "published_age_length_ma": 0.0, "tolerance_ma": tolerance})
            node.dist = 0.0
    return (newick_text(tree) if adjustments else text), adjustments


def date(tree, provenance, calibrations, outdir, command, settings=None, threads=1):
    cfg = validate_settings({} if settings is None else settings)
    if type(threads) is not int or threads != 1:
        raise ValueError("the LSD2 workflow requires threads=1")
    source_tree, tree = tree, read_tree(tree)
    source_qc = json.loads(Path(provenance).read_text())
    if source_qc.get("branch_length_unit") != "substitutions_per_site" or not source_qc.get("outgroup"):
        raise ValueError("dating requires substitution-per-site lengths and a selected basal group")
    outgroup = outgroup_ids(source_qc["outgroup"])
    # Orient a copy for MRCA calibration checks; LSD2 itself roots its input below.
    root_on_outgroup(tree, outgroup)
    numsites = source_qc.get("total_gene_sites")
    if type(numsites) is not int or not 0 < numsites < 2**31:
        raise ValueError("LSD2 requires total_gene_sites as a positive 32-bit integer in species-tree provenance")
    if not any(n.dist > 0 for n in tree.traverse() if not n.is_root):
        raise ValueError("dating requires positive substitution lengths")
    if any(not re.fullmatch(r"[A-Za-z0-9_.-]+", n) for n in tree.leaf_names()):
        raise ValueError("LSD2 tip labels must contain only letters, digits, underscores, dots or hyphens")
    used = set(tree.leaf_names())
    for number, node in enumerate(tree.traverse(), 1):
        if not node.is_leaf:
            node.name = f"N{number:06d}"
            while node.name in used:
                node.name = "N" + node.name
            used.add(node.name)
    constraints = calibration_rows(tree, calibrations)
    command = executable(command)
    build, binary = lsd2_build(command)
    out = Path(outdir).resolve()
    runs = out / "lsd2_runs"
    runs.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix="run-", dir=runs))
    # Retain the supplied tree. MRCA constraints are resolved after LSD2 rooting.
    (work / "lsd2.input.nwk").write_text(Path(source_tree).read_text())
    (work / "lsd2.dates.txt").write_text(dates_text(constraints))
    (work / "lsd2.outgroups.txt").write_text(str(len(outgroup)) + "\n" + "\n".join(outgroup) + "\n")
    argv = [command, "-i", "lsd2.input.nwk", "-d", "lsd2.dates.txt", "-g", "lsd2.outgroups.txt",
            "-o", "lsd2.result", "-r", "k", "-z", "0", "-s", str(numsites), "-l", "-1",
            "-u", "0", "-U", "0", "-D", "1", "-v", str(cfg["variance"])]
    if cfg["variance_parameter"] is not None:
        argv += ["-b", str(cfg["variance_parameter"])]
    invocation = {"argv": argv, "cwd": str(work)}
    write_json(work / "lsd2.command.json", invocation)
    print(f"LSD2 dating on the selected basal edge: {work}", flush=True)
    started = time.monotonic()
    with (work / "lsd2.log").open("w") as log:
        subprocess.run(argv, cwd=work, stdout=log, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL,
                       check=True, env={**os.environ, "OMP_NUM_THREADS": "1"})
    elapsed = time.monotonic() - started
    result = native_result(work / "lsd2.result", cfg["variance"])
    native_text = extract_time_tree(work / "lsd2.result.date.nexus")
    (work / "lsd2.dated.nwk").write_text(native_text)
    checked_text, adjustments = normalize_numerical_zeros(native_text, max(result["root_age_solution_bounds_ma"]))
    (work / "validation.nwk").write_text(checked_text)
    dated, ages, height = validate_dated(work / "validation.nwk", tree, constraints)
    validate_root(dated, outgroup)
    # Publish only validated native time edges; restore stable pipeline node names.
    for path in work.iterdir():
        if path.is_file() and path.name.startswith("lsd2."):
            with atomic_writer(out / path.name) as handle:
                handle.write(path.read_text())
    with atomic_writer(out / "species_tree.dated.nwk") as handle:
        handle.write(newick_text(dated))
    write_tsv(out / "calibrations.resolved.tsv", list(constraints[0]), constraints)
    write_tsv(out / "node_ages.tsv", ["node", "age_ma"],
              [{"node": node, "age_ma": age} for node, age in ages.items()])
    write_tsv(out / "numerical_zero_adjustments.tsv",
              ["preorder_index", "native_age_length_ma", "published_age_length_ma", "tolerance_ma"], adjustments)
    write_json(out / "provenance.json", {
        "created_at": now(), "commands": [invocation], "executable": binary, "build": build,
        "tree": file_record(source_tree), "source_provenance": file_record(provenance),
        "calibrations": file_record(calibrations), "branch_length_unit": "million_years",
        "method": "LSD2 least squares", "rate_model": "least squares with a shared rate and branch residuals",
        "settings": cfg, "numsites": numsites, "numsites_source": "sum of retained trimmed gene sites",
        "numsites_is_effective_sample_size": False, "confidence_intervals": False,
        "outgroup": outgroup, "outgroups_retained": True, "root_edge_preserved": True,
        "root_position_reestimated": True, "topology_preserved": True, "tip_age_ma": 0,
        "short_branches_collapsed": False, "outlier_tips_removed": False,
        "concatenation_used": False, "root_age_ma": height, "workdir": str(work),
        "elapsed_seconds": elapsed, "native_time_branch_lengths_preserved": not adjustments,
        "numerical_zero_adjustments": len(adjustments),
        "node_age_calculation": "mean descendant-tip distance from native rounded time lengths",
        "point_date_policy": "native estimate; native midpoint of boundary solutions when time scale is nonunique",
        "solution_bounds_are_confidence_intervals": False,
        "optimization_validation": "native completion, finite solution, topology, temporal and calibration checks",
        **result, "review_status": "needs_review" if result["diagnostics"] or not result["unique_time_scale"] else "checks_passed"})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("tree", "provenance", "calibrations", "outdir", "command"):
        parser.add_argument("--" + name, required=True)
    parser.add_argument("--settings", type=json.loads, default={})
    parser.add_argument("--threads", type=int, default=1)
    date(**vars(parser.parse_args()))
