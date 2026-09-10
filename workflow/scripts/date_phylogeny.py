#!/usr/bin/env python3
"""Date CASTLES substitution lengths with unmodified, standalone LSD2."""
import argparse
import json
import math
import os
from pathlib import Path
import re
import statistics
import subprocess
import tempfile

from common import atomic_writer, file_record, now, read_tsv, write_json, write_tsv
from infer_phylogeny import executable, read_tree

DEFAULT_SETTINGS = {"variance": 1, "variance_parameter": None, "numsites": None}


def validate_settings(settings):
    if not isinstance(settings, dict):
        raise ValueError("LSD2 settings must be a mapping")
    unknown = set(settings) - set(DEFAULT_SETTINGS)
    if unknown:
        raise ValueError(f"unknown LSD2 settings: {sorted(unknown)}")
    cfg = {**DEFAULT_SETTINGS, **settings}
    if type(cfg["variance"]) is not int or cfg["variance"] not in {0, 1, 2}:
        raise ValueError("LSD2 variance must be 0, 1, or 2")
    value = cfg["variance_parameter"]
    if value is not None and (type(value) not in {int, float} or not math.isfinite(value) or value <= 0):
        raise ValueError("LSD2 variance_parameter must be null or finite and positive")
    if cfg["variance"] == 0 and value is not None:
        raise ValueError("LSD2 variance_parameter requires variance 1 or 2")
    if cfg["numsites"] is not None and (type(cfg["numsites"]) is not int or not 0 < cfg["numsites"] < 2**31):
        raise ValueError("LSD2 numsites must be null or a positive 32-bit integer")
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


def validate_dated(path, original, constraints):
    dated = read_tree(path, original.leaf_names())
    # LSD2 may omit/replace internal labels. Match children by descendant
    # signatures, in linear memory, and restore our stable node names.
    # Intern signatures jointly to avoid quadratic storage on ladder-like trees.
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
    names = {code: n.name for n, code in a.items()}
    for node, code in b.items():
        node.name = names[code]
    depths = {dated: 0.0}
    for node in dated.traverse():
        if not node.is_root:
            depths[node] = depths[node.up] + node.dist
    tips = [depths[n] for n in dated.leaves()]
    height = statistics.mean(tips)
    tolerance = max(1e-8, height * 1e-7)
    if not math.isfinite(height) or height <= 0 or max(tips) - min(tips) > tolerance:
        raise ValueError("LSD2 time tree is not ultrametric for extant tips")
    ages = {n.name: height - d for n, d in depths.items() if not n.is_leaf}
    for row in constraints:
        age = ages[row["node"]]
        if not row["min_age_ma"] - tolerance <= age <= row["max_age_ma"] + tolerance:
            raise ValueError(f"LSD2 output violates calibration at {row['node']}")
    return dated, ages, height


def newick_text(tree):
    from ete4.parser.newick import make_parser
    return tree.write(parser=make_parser(1, dist="%.17g"), format_root_node=True) + "\n"


def lsd_dates(tree, constraints):
    """LSD2 dates increase toward the present; ages in Ma increase backward."""
    lines = [f"{name} 0" for name in sorted(tree.leaf_names())]
    for row in constraints:
        low, high = -row["max_age_ma"], -row["min_age_ma"]
        value = f"{low:.17g}" if low == high else f"b({low:.17g},{high:.17g})"
        # MRCA expressions also work for the root, whose label upstream may omit.
        lines.append(f"mrca({row['taxa']}) {value}")
    return str(len(lines)) + "\n" + "\n".join(lines) + "\n"


def parse_lsd_report(path, variance=1):
    text = Path(path).read_text()
    rows = re.findall(r"^\s*rate\s+([^,\s]+)\s*,\s*tMRCA\s+([^,\s]+)\s*,\s*objective function\s+(\S+)\s*$",
                      text, flags=re.MULTILINE)
    expected = 2 if variance == 2 else 1
    if len(rows) != expected or (variance == 2 and "Results of the second run" not in text):
        raise ValueError("LSD2 did not report the expected completed single-rate fit(s)")
    def interval(token):
        values = list(map(float, token.split(":")))
        if len(values) == 1:
            values *= 2
        if len(values) != 2 or not all(math.isfinite(v) for v in values) or values[0] > values[1]:
            raise ValueError("LSD2 reported an invalid/nonfinite solution interval")
        return values
    fits = []
    for rate_token, date_token, objective_token in rows:
        rate, root_date, objective = interval(rate_token), interval(date_token), float(objective_token)
        if rate[0] <= 0 or root_date[1] >= 0 or not math.isfinite(objective) or objective < 0:
            raise ValueError("LSD2 reported an invalid/nonfinite rate, root date or objective")
        fits.append({"rate_interval": rate, "root_date_interval": root_date, "objective": objective})
    version = re.search(r"LEAST-SQUARE.*?v\.([0-9.]+)", text)
    if not version or version[1] != "2.4.4":
        raise ValueError("unsupported LSD2 report version; use workflow/envs/dating.yaml")
    warning = re.search(r"\*WARNINGS:\s*(.*?)\*RESULTS:", text, flags=re.DOTALL)
    final = fits[-1]
    return {"version": version[1], **final, "fits": fits,
            "rate_substitutions_per_site_per_ma": final["rate_interval"][0] if len(set(final["rate_interval"])) == 1 else None,
            "unique_scale": len(set(final["rate_interval"])) == 1,
            "warnings": warning[1].strip().splitlines() if warning else []}


def rounding_tolerance(value):
    # Upstream LSD2 2.4.4 uses default C++ stream precision: six significant digits.
    # Include omitted trailing zeroes, e.g. "100" represents rounding at 0.001.
    return max(1e-12, 0.500001 * 10 ** (math.floor(math.log10(abs(value))) - 5)) if value else 1e-12


def read_lsd_dates(path, original, constraints):
    """Validate native dates, then export time branches without accumulated rounding."""
    from Bio import Phylo
    native = Phylo.read(path, "nexus")
    tips = [n.name for n in native.get_terminals()]
    if len(tips) != len(set(tips)) or set(tips) != set(original.leaf_names()):
        raise ValueError("LSD2 changed the species set")
    # Intern descendant signatures jointly. Storage remains linear on ladder trees.
    intern, source_codes, native_codes = {}, {}, {}
    for node in original.traverse("postorder"):
        key = ("tip", node.name) if node.is_leaf else ("node", tuple(sorted(source_codes[c] for c in node.children)))
        source_codes[node] = intern.setdefault(key, len(intern))
    for node in native.find_clades(order="postorder"):
        key = ("tip", node.name) if node.is_terminal() else ("node", tuple(sorted(native_codes[c] for c in node.clades)))
        native_codes[node] = intern.setdefault(key, len(intern))
    if source_codes[original] != native_codes[native.root] or len(source_codes) != len(native_codes):
        raise ValueError("LSD2 changed the rooted topology")
    names = {code: node.name for node, code in source_codes.items()}
    raw, tolerances = {}, {}
    for node, code in native_codes.items():
        matches = re.findall(r'(?:^|[&,\s])date=(?:"([^"]+)"|([^,\]\s]+))', node.comment or "")
        if len(matches) != 1:
            raise ValueError("LSD2 native time tree has missing/duplicate node dates")
        age = -float(matches[0][0] or matches[0][1])
        if not math.isfinite(age) or age < 0:
            raise ValueError("LSD2 returned a nonfinite or future node age")
        if node.is_terminal() and age != 0:
            raise ValueError("LSD2 tip dates must all be zero for extant species")
        name = names[code]
        raw[name], tolerances[name] = age, rounding_tolerance(age)
    for node in native.find_clades():
        parent = names[native_codes[node]]
        for child in node.clades:
            name = names[native_codes[child]]
            length = child.branch_length
            if length is None or not math.isfinite(length) or length < 0 or raw[parent] < raw[name]:
                raise ValueError("LSD2 violated temporal constraints or returned invalid time branches")
            tolerance = tolerances[parent] + tolerances[name] + rounding_tolerance(length)
            if abs(length - (raw[parent] - raw[name])) > tolerance:
                raise ValueError("LSD2 native branches do not match node dates in time units")
    bounds = {row["node"]: row for row in constraints}
    lower, upper = {}, {}
    for node in original.traverse("postorder"):
        name = node.name
        lo, hi = (0.0, 0.0) if node.is_leaf else (max(0, raw[name] - tolerances[name]), raw[name] + tolerances[name])
        if name in bounds:
            lo = max(lo, bounds[name]["min_age_ma"])
            hi = min(hi, bounds[name]["max_age_ma"])
        lo = max(lo, max((lower[c.name] for c in node.children), default=0))
        if lo > hi:
            raise ValueError(f"LSD2 output violates calibration/temporal bounds beyond rounding at {name}")
        lower[name], upper[name] = lo, hi
    dated = original.copy()
    ages, adjustments = {}, []
    for node in dated.traverse():
        name = node.name
        hi = min(upper[name], ages[node.up.name]) if not node.is_root else upper[name]
        age = min(max(raw[name], lower[name]), hi)
        ages[name] = age
        if age != raw[name]:
            adjustments.append({"node": name, "native_age_ma": raw[name], "age_ma": age,
                                "adjustment_ma": age - raw[name], "rounding_tolerance_ma": tolerances[name]})
        node.dist = ages[node.up.name] - age if not node.is_root else 0
    return dated, ages, raw, adjustments


def date(tree, provenance, calibrations, outdir, command, settings=None, threads=1):
    cfg = validate_settings({} if settings is None else settings)
    if type(threads) is not int or threads != 1:
        raise ValueError("standalone LSD2 requires threads=1")
    source_tree, tree = tree, read_tree(tree)
    source_qc = json.loads(Path(provenance).read_text())
    if source_qc.get("branch_length_unit") != "substitutions_per_site" or not source_qc.get("outgroup"):
        raise ValueError("dating requires a rooted tree with substitution-per-site branch lengths")
    if source_qc["outgroup"] not in set(tree.leaf_names()):
        raise ValueError("dating outgroup is absent from the species tree")
    numsites = cfg["numsites"] if cfg["numsites"] is not None else source_qc.get("total_gene_sites")
    if type(numsites) is not int or not 0 < numsites < 2**31:
        raise ValueError("LSD2 requires total_gene_sites in species-tree provenance (or explicit lsd2.numsites)")
    if len(tree.children) != 2:
        raise ValueError("LSD2 requires an explicitly rooted species tree")
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
    out = Path(outdir).resolve()
    runs = out / "lsd2_runs"
    runs.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix="run-", dir=runs))
    input_text = newick_text(tree)
    (work / "input.nwk").write_text(input_text)
    (work / "dates.txt").write_text(lsd_dates(tree, constraints))
    argv = [command, "-i", "input.nwk", "-d", "dates.txt", "-o", "dated",
            "-s", str(numsites), "-l", "-1", "-u", "0", "-U", "0", "-v", str(cfg["variance"])]
    if cfg["variance_parameter"] is not None:
        argv += ["-b", f"{cfg['variance_parameter']:.17g}"]
    invocation = {"argv": argv, "cwd": str(work)}
    write_json(work / "command.json", invocation)
    print(f"LSD2: {work}", flush=True)
    with (work / "run.log").open("w") as log:
        subprocess.run(argv, cwd=work, stdout=log, stderr=subprocess.STDOUT, check=True,
                       env={**os.environ, "OMP_NUM_THREADS": "1"})
    report = parse_lsd_report(work / "dated", cfg["variance"])
    dated, ages, raw, adjustments = read_lsd_dates(work / "dated.date.nexus", tree, constraints)
    root_low, root_high = report["root_date_interval"]
    root_tolerance = 2 * max(rounding_tolerance(root_low), rounding_tolerance(root_high))
    if not root_low - root_tolerance <= -ages[tree.name] <= root_high + root_tolerance:
        raise ValueError("LSD2 report and time tree have inconsistent root dates")
    normalized = work / "species_tree.dated.nwk"
    normalized.write_text(newick_text(dated))
    validate_dated(normalized, tree, constraints)
    # Publish only after validation. The final provenance file is the completion record.
    files = {"species_tree.dated.nwk": normalized, "lsd2.input.nwk": work / "input.nwk",
             "lsd2.dates.txt": work / "dates.txt", "lsd2.report.txt": work / "dated",
             "lsd2.dated.date.nexus": work / "dated.date.nexus", "lsd2.fitted.nwk": work / "dated.nwk",
             "lsd2.command.json": work / "command.json"}
    if any(not p.is_file() or p.stat().st_size == 0 for p in files.values()):
        raise ValueError("LSD2 did not produce all required outputs")
    build_path = Path(command).parent.parent / "share" / "lsd2" / "build.json"
    build = json.loads(build_path.read_text()) if build_path.is_file() else None
    executable_record = file_record(command)
    if build and build.get("executable_sha256") != executable_record["sha256"]:
        raise ValueError("LSD2 executable does not match its build provenance")
    for name, path in files.items():
        with atomic_writer(out / name, "wb") as handle:
            handle.write(path.read_bytes())
    rows = [{"node": n.name, "age_ma": ages[n.name], "native_age_ma": raw[n.name],
             "rounding_adjustment_ma": ages[n.name] - raw[n.name]} for n in dated.traverse() if not n.is_leaf]
    write_tsv(out / "node_ages.tsv", list(rows[0]), rows)
    write_tsv(out / "calibrations.resolved.tsv", list(constraints[0]), constraints)
    write_tsv(out / "rounding_adjustments.tsv",
              ["node", "native_age_ma", "age_ma", "adjustment_ma", "rounding_tolerance_ma"], adjustments)
    diagnostics = list(report["warnings"])
    if report["rate_interval"][0] <= 1e-10:
        diagnostics.append("estimated rate is at LSD2's lower bound")
    write_json(out / "provenance.json", {
        "created_at": now(), "commands": [invocation], "executable": executable_record, "build": build,
        "tree": file_record(source_tree), "source_provenance": file_record(provenance),
        "calibrations": file_record(calibrations), "branch_length_unit": "million_years",
        "method": "LSD2 least-squares dating", "rate_model": "single estimated substitution rate",
        "settings": cfg, "numsites": numsites,
        "numsites_source": "explicit override" if cfg["numsites"] is not None else "sum of retained trimmed gene sites",
        "numsites_is_effective_sample_size": False, "confidence_intervals": False,
        "root_preserved": True, "topology_preserved": True, "tip_age_ma": 0,
        "native_report": report, "root_age_ma": ages[tree.name], "concatenation_used": False,
        "scale_selection": "unique native scale" if report["unique_scale"] else "native midpoint of feasible date bounds; not confidence intervals",
        "native_date_significant_digits": 6,
        "time_branch_export": "parent age minus child age, written with 17 significant digits",
        "rounding_adjustments": len(adjustments),
        "rounding_policy": "intersect native rounding intervals with hard calibrations and ancestor age bounds",
        "diagnostics": diagnostics, "review_status": "needs_review" if diagnostics else "checks_passed",
        "native_workdir": str(work)})
    for message in diagnostics:
        print("LSD2 diagnostic: " + message, flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["tree", "provenance", "calibrations", "outdir", "command"]:
        parser.add_argument("--" + name, required=True)
    parser.add_argument("--settings", default="{}")
    parser.add_argument("--threads", type=int, default=1)
    args = vars(parser.parse_args())
    args["settings"] = json.loads(args["settings"])
    date(**args)
