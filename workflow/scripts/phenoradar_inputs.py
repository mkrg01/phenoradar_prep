#!/usr/bin/env python3
"""Collect completed outputs for present and future PhenoRadar analyses."""
import argparse
import csv
import gzip
import json
import math
from pathlib import Path
import tempfile

from common import read_tsv, sha256, species_from_gene_id
from export_species_tpm import export as export_species_tpm
from filter_species import fasta_records, manifest, safe_name, table, validate_exclusions
from phenoradar_metadata import read_base, with_pairs
from layout import (ORTHOGROUP_MAPPING, ORTHOGROUP_EXPRESSION, ORTHOGROUP_ALIGNMENTS,
                    PHYLOGENY_BRANCHES, REPRESENTATIVES)


def discover(source, exclusions=()):
    """Freeze existing files as absolute inputs, without requesting producers."""
    exclusions = validate_exclusions(list(exclusions))
    original = Path(source).resolve()
    _, original_samples, original_species, _ = manifest(original / "metadata/samples.tsv")
    if set(exclusions) - set(original_species):
        raise ValueError("exclude_species contains species outside the original samples")
    source, filtered = original, None
    files = {original / "metadata/samples.tsv"}
    if exclusions:
        source = original / "filtered"
        marker = source / "manifest.json"
        if not marker.is_file():
            raise ValueError("run filter_species with the current exclude_species before phenoradar_inputs")
        filtered = json.loads(marker.read_text())
        retained = set(original_species) - set(exclusions)
        retained_runs = {r["run"] for r in original_samples if r["species"] in retained}
        if (filtered.get("report_type") != "species_filter" or filtered.get("source") != str(original)
                or filtered.get("exclude_species") != exclusions
                or set(filtered.get("retained_species", [])) != retained
                or set(filtered.get("retained_runs", [])) != retained_runs or not retained):
            raise ValueError("filtered snapshot does not match the current source/exclude_species; rerun filter_species")
        files.add(marker)
    for name in ["metadata/samples.tsv", "metadata/species_metadata.tsv"]:
        path = source / name
        if not path.is_file():
            raise ValueError(f"missing required workflow output: {path}")
        files.add(path)

    sections, links, hashes, trees, contrasts = {}, {}, {}, {}, []

    def publish(relative, destination=None):
        path = source / relative
        if path.is_file():
            files.add(path)
            links[destination or relative] = path

    def folder(relative, recursive=False):
        root = source / relative
        for path in sorted(root.rglob("*") if recursive else root.glob("*")):
            if path.is_file() and not any(p.startswith(".") for p in path.relative_to(root).parts):
                publish(str(path.relative_to(source)))

    def branch(name, required):
        present = [p for p in required if (source / p).is_file()]
        status = "ready" if len(present) == len(required) else "incomplete" if present else "absent"
        sections[name] = status
        files.update(source / p for p in present)
        return status == "ready"

    folder("metadata")
    publish("run.json")
    publish("metadata/species_metadata.tsv", "species_metadata.tsv")
    if branch("orthogroups", [f"{ORTHOGROUP_EXPRESSION}/tpm.tsv", f"{ORTHOGROUP_EXPRESSION}/mapping_qc.tsv"]):
        folder(ORTHOGROUP_EXPRESSION)
        publish(f"{ORTHOGROUP_EXPRESSION}/tpm.tsv", "tpm.tsv")
    if branch("mapping", [f"{ORTHOGROUP_MAPPING}/{name}" for name in
                          ["gene_orthogroups.tsv", "mappings.sqlite", "merge_qc.json"]]):
        folder(ORTHOGROUP_MAPPING)
    if branch("kegg", ["kegg/ko_tpm_sum.tsv", "kegg/mapping_qc.tsv", "kegg/ko_support.tsv"]):
        for name in ["ko_tpm_sum.tsv", "ko_tpm_sum_wide.tsv", "ko_support.tsv", "mapping_qc.tsv",
                     "genes.tsv", "gene_kos.tsv"]:
            publish("kegg/" + name)
    for group in ["module", "pathway"]:
        relative = f"kegg/ko_{group}s.tsv"
        if branch(f"kegg_{group}s", [relative, "kegg/reference_qc.json"]):
            publish(relative)
            publish("kegg/reference_qc.json")
            record = json.loads((source / "kegg/reference_qc.json").read_text())
            maps = {Path(r["path"]).name: r["sha256"] for r in record.get("results", [])}
            if Path(relative).name not in maps:
                raise ValueError("KO membership is not recorded in KEGG reference completion: " + relative)
            hashes[str(source / relative)] = maps[Path(relative).name]

    completion = f"{ORTHOGROUP_ALIGNMENTS}/filter_qc.json" if filtered else f"{ORTHOGROUP_ALIGNMENTS}/provenance.json"
    if branch("alignments", [completion]):
        report = json.loads((source / completion).read_text())
        seen = set()
        for record in report["alignments"]:
            og = safe_name(record["orthogroup"])
            if og in seen:
                raise ValueError("duplicate OG in alignment completion record: " + og)
            seen.add(og)
            if filtered and record["after"] == 0:
                continue
            relative = f"{ORTHOGROUP_ALIGNMENTS}/{og}.faa"
            path = source / relative
            if not path.is_file():
                raise ValueError("completed alignment is missing: " + str(path))
            publish(relative)
            publish(relative, f"alignments/{og}.faa")
            if not filtered:
                hashes[str(path)] = record["alignment"]["sha256"]
        for path in (source / ORTHOGROUP_ALIGNMENTS).glob("*"):
            if path.suffix != ".faa":
                publish(str(path.relative_to(source)))
    elif any((source / ORTHOGROUP_ALIGNMENTS).glob("*.faa")):
        sections["alignments"] = "incomplete"

    # Species/run identity comes from the selected manifest, never stale files.
    _, samples, _, _ = manifest(source / "metadata/samples.tsv")
    for name in sorted({r["odb_species"] for r in samples}):
        protein = f"proteins/{name}_protein.fa"
        report = f"proteins/{name}_protein.json"
        if branch(f"proteins/{name}", [protein] if filtered else [protein, report]):
            publish(protein)
            publish(report)
            if not filtered:
                hashes[str(source / protein)] = json.loads((source / report).read_text())["protein"]["sha256"]
        relative = f"kegg/species/{name}"
        if branch(relative, [f"{relative}/{n}" for n in ["genes.tsv", "gene_kos.tsv", "detail.tsv", "provenance.json"]]):
            folder(relative)

    for relative in [*PHYLOGENY_BRANCHES.values(), REPRESENTATIVES]:
        if filtered and branch(relative, [f"{relative}/pruning.json"]):
            folder(relative, recursive=True)
            tree = source / relative / "species_tree.pruned.nwk"
            if tree.is_file():
                trees[str(tree)] = False
            dated = source / relative / "dating/species_tree.dated.pruned.nwk"
            if dated.is_file():
                trees[str(dated)] = False
        else:
            for name, marker in [("species_tree.nwk", "species_tree.json"),
                                 ("gene_trees.nwk", "gene_trees.json")]:
                if branch(f"{relative}/{name}", [f"{relative}/{name}", f"{relative}/{marker}"]):
                    publish(f"{relative}/{name}")
                    publish(f"{relative}/{marker}")
                    publish(f"{relative}/species_coverage.tsv")
                    if name == "species_tree.nwk":
                        trees[str(source / relative / name)] = relative != PHYLOGENY_BRANCHES["all"]
                    else:
                        record = json.loads((source / relative / marker).read_text())
                        for locus in record.get("retained", []):
                            og = safe_name(locus["marker"])
                            for directory in ["markers", "alignments", "alignments/raw", "gene_trees"]:
                                for suffix in [".faa", ".nwk", ".json", ".columns.tsv"]:
                                    publish(f"{relative}/{directory}/{og}{suffix}")
            for directory, marker in [("selection", "selection.json"), ("plan", "provenance.json"),
                                      ("rooting", "outgroup.json"), ("timetree", "provenance.json"),
                                      ("dating", "provenance.json"), ("taxonomy_check", "summary.json")]:
                required = [f"{relative}/{directory}/{marker}"]
                if directory == "dating":
                    required.append(f"{relative}/dating/species_tree.dated.nwk")
                if branch(f"{relative}/{directory}", required):
                    folder(f"{relative}/{directory}", recursive=directory == "taxonomy_check")
                    if directory == "dating":
                        trees[str(source / relative / "dating/species_tree.dated.nwk")] = relative != PHYLOGENY_BRANCHES["all"]
            # Completed sequence/alignment steps are useful even before tree inference.
            if (source / relative / "plan/provenance.json").is_file():
                for directory, field, plan, suffix in [
                    ("species", "species", "species", ".faa"),
                    ("alignments/raw", "marker", "markers", ".faa"),
                    ("alignments", "marker", "markers", ".faa"),
                    ("gene_trees", "marker", "markers", ".nwk"),
                ]:
                    plan_path = source / relative / "plan" / f"{plan}.tsv"
                    if not plan_path.is_file():
                        continue
                    for row in read_tsv(plan_path):
                        stem = f"{relative}/{directory}/{safe_name(row[field])}"
                        if all((source / (stem + ext)).is_file() for ext in [suffix, ".json"]):
                            for ext in [suffix, ".json", ".columns.tsv"]:
                                publish(stem + ext)
                if (source / relative / "markers/.snakemake_timestamp").is_file():
                    files.add(source / relative / "markers/.snakemake_timestamp")
                    plan_path = source / relative / "plan/markers.tsv"
                    if plan_path.is_file():
                        for row in read_tsv(plan_path):
                            publish(f"{relative}/markers/{safe_name(row['marker'])}.faa")
        contrast = f"{relative}/contrast"
        if branch(contrast, [f"{contrast}/species_metadata.tsv", f"{contrast}/summary.json"]):
            folder(contrast)
            contrasts.append(source / contrast / "species_metadata.tsv")

    for directory in ["", "orthogroups", ORTHOGROUP_MAPPING]:
        for suffix in [".tsv", ".tsv.gz"]:
            publish(str(Path(directory) / ("orthogroup_annotations" + suffix)))
    return dict(source=str(source), files=sorted(str(p) for p in files), links=links,
                sections=sections, filtered=filtered, contrasts=contrasts, hashes=hashes, trees=trees,
                original_runs={r["run"]: r["species"] for r in original_samples if r["species"] not in exclusions})


def validate_expression(path, runs, feature, value):
    """Validate run-level values without choosing or combining biological replicates."""
    reader = table(path, ["species", "run", feature, value])
    next(reader)
    seen, features, current_features = set(), set(), set()
    current = None
    for row in reader:
        run = row["run"]
        if runs.get(run) != row["species"]:
            raise ValueError(f"expression run/species differs from samples: {path}: {row['species']}/{run}")
        if run != current:
            if run in seen:
                raise ValueError(f"duplicate/discontiguous run block in expression: {path}: {run}")
            seen.add(run)
            current, current_features = run, set()
        name = safe_name(row[feature])
        if name in current_features:
            raise ValueError(f"duplicate run/feature coordinate: {path}: {run}/{name}")
        current_features.add(name)
        features.add(name)
        try:
            number = float(row[value])
        except ValueError as error:
            raise ValueError(f"invalid numeric expression: {path}: {row[value]!r}") from error
        if not math.isfinite(number) or number < 0:
            raise ValueError(f"expression must be finite and nonnegative: {path}")
    if seen != set(runs):
        raise ValueError(f"expression lacks selected runs: {path}: {sorted(set(runs) - seen)}")
    return features


def validate_qc(path, runs):
    reader = table(path, ["species", "run"])
    next(reader)
    seen = set()
    for row in reader:
        run = row["run"]
        if run in seen or runs.get(run) != row["species"]:
            raise ValueError("expression QC run/species differs from samples: " + str(path))
        seen.add(run)
    if seen != set(runs):
        raise ValueError("expression QC lacks selected runs: " + str(path))


def validate_groups(path, group):
    reader = table(path, ["ko", group])
    next(reader)
    seen = set()
    for row in reader:
        pair = (safe_name(row["ko"]), safe_name(row[group]))
        if pair in seen:
            raise ValueError("duplicate KO/group membership: " + str(path))
        seen.add(pair)


def validate_alignment(path, species):
    seen, width = set(), None
    for header, sequence in fasta_records(path):
        gene = header.split()[0]
        if species_from_gene_id(gene) not in species or gene in seen:
            raise ValueError("unexpected species or duplicate gene in alignment: " + str(path))
        seen.add(gene)
        if width is not None and len(sequence) != width:
            raise ValueError("unequal alignment lengths: " + str(path))
        width = len(sequence)
    if not seen:
        raise ValueError("empty alignment: " + str(path))


def validate_tree(path, species, subset=False):
    from ete4 import Tree
    text = path.read_text().strip()
    if text.count(";") != 1:
        raise ValueError("expected one Newick species tree")
    tree = Tree(text, parser=1)
    names = list(tree.leaf_names())
    if (len(names) != len(set(names)) or not names or not set(names) <= species
            or (not subset and set(names) != species)):
        raise ValueError("tree species differ from selected metadata; use the matching full/pruned tree")
    for node in tree.traverse():
        if node.dist is not None and (not math.isfinite(node.dist) or node.dist < 0):
            raise ValueError("tree branch lengths must be finite and nonnegative when present")


def validate_annotations(path, orthogroups):
    matched = set()
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) != 3 or not row[0] or not row[1].isdigit():
                raise ValueError("OG annotations require headerless OG/taxid/description rows: " + str(path))
            if row[0] in orthogroups:
                if row[0] in matched:
                    raise ValueError("duplicate annotation for selected OG: " + row[0])
                matched.add(row[0])
    if orthogroups and not matched:
        raise ValueError("OG annotations do not match expression IDs; check the OrthoDB reference version")
    print(f"OG descriptions: {len(matched)}/{len(orthogroups)} expressed OGs", flush=True)


def check_destination(out):
    if out.is_symlink():
        raise ValueError("phenoradar_inputs output directory must not be a symlink")
    allowed = {"species_metadata.tsv", "tpm.tsv", "orthogroup_annotations.tsv",
               "orthogroup_annotations.tsv.gz", "run.json", ".snakemake_timestamp"}
    directories = {"metadata", "proteins", "orthogroups", "phylogeny", "alignments", "kegg"}
    if out.exists():
        for path in out.iterdir():
            if path.name in allowed and (path.is_symlink() or
                    (path.name in {"species_metadata.tsv", "tpm.tsv", ".snakemake_timestamp"} and path.is_file())):
                continue
            if path.name in directories and path.is_dir() and not path.is_symlink():
                for child in path.rglob("*"):
                    if child.is_symlink() and not child.is_dir():
                        continue
                    if child.is_dir() and not child.is_symlink():
                        continue
                    raise ValueError("unrecognized file in output directory: " + str(child))
                continue
            raise ValueError("refusing to replace unrecognized output contents: " + str(path))


def export(source, outdir=None, exclusions=(), trait="C4"):
    out = Path(outdir) if outdir else Path(source) / "phenoradar_inputs"
    check_destination(out)
    out = out.resolve()
    inventory = discover(source, exclusions)
    source = Path(inventory["source"])
    paths = [Path(p) for p in inventory["files"]]
    if out == source or any(out == p.resolve() or out in p.resolve().parents for p in paths):
        raise ValueError("phenoradar_inputs output overlaps its source files")
    stats = {p: (p.stat().st_size, p.stat().st_mtime_ns, p.stat().st_ino) for p in paths}
    print("PhenoRadar source: " + str(source), flush=True)
    for name, status in inventory["sections"].items():
        print(f"{name}: {status}", flush=True)
    expected_hashes = dict(inventory["hashes"])
    filtered = inventory["filtered"]
    if filtered:
        recorded = {r["path"]: r["sha256"] for r in filtered["outputs"]}
        for path in paths:
            if path.is_relative_to(source) and path != source / "manifest.json":
                relative = str(path.relative_to(source))
                if relative not in recorded:
                    raise ValueError("input is not recorded in completed filtered snapshot: " + relative)
                if str(path) in expected_hashes and expected_hashes[str(path)] != recorded[relative]:
                    raise ValueError("filtered checksum differs from original completion: " + relative)
                expected_hashes[str(path)] = recorded[relative]
    for path, expected in expected_hashes.items():
        if sha256(path) != expected:
            raise ValueError("completed input checksum mismatch: " + path)
    _, samples, by_species, runs = manifest(source / "metadata/samples.tsv")
    species = set(by_species)
    if filtered and species != set(filtered["retained_species"]):
        raise ValueError("filtered samples differ from recorded retained species")
    if filtered and {r["run"]: r["species"] for r in samples} != inventory["original_runs"]:
        raise ValueError("filtered runs differ from the current original samples; rerun filter_species")
    metadata = read_base(source / "metadata/species_metadata.tsv", trait)
    if {r["species"] for r in metadata} != species:
        raise ValueError("PhenoRadar metadata species differ from samples")
    for pairs in inventory["contrasts"]:
        with_pairs(metadata, pairs, trait)
    links, orthogroups = inventory["links"], set()
    for dest, feature, value, qc in [("tpm.tsv", "orthogroup", "tpm", f"{ORTHOGROUP_EXPRESSION}/mapping_qc.tsv"),
                                    ("kegg/ko_tpm_sum.tsv", "ko", "tpm_sum", "kegg/mapping_qc.tsv")]:
        if dest in links:
            print("Checking expression: " + str(links[dest]), flush=True)
            features = validate_expression(links[dest], runs, feature, value)
            validate_qc(source / qc, runs)
            if feature == "orthogroup":
                orthogroups = features
    for dest, path in links.items():
        if dest.startswith("alignments/"):
            validate_alignment(path, species)
        elif dest in {"kegg/ko_modules.tsv", "kegg/ko_pathways.tsv"}:
            validate_groups(path, "module" if "modules" in dest else "pathway")
        elif str(path) in inventory["trees"]:
            validate_tree(path, species, subset=inventory["trees"][str(path)])
        elif path.name.startswith("orthogroup_annotations."):
            validate_annotations(path, orthogroups)
    out.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".phenoradar-inputs-", dir=out.parent) as tmp:
        stage = Path(tmp) / "inputs"
        stage.mkdir()
        for dest, path in links.items():
            target = stage / dest
            target.parent.mkdir(parents=True, exist_ok=True)
            if dest == "tpm.tsv":
                export_species_tpm(source / "metadata/samples.tsv", path, target)
                print(f"{dest}: species/orthogroup/tpm exported from {path}", flush=True)
            else:
                target.symlink_to(path.resolve())
                print(f"{dest} -> {path}", flush=True)
        for path, expected in stats.items():
            stat = path.stat()
            if (stat.st_size, stat.st_mtime_ns, stat.st_ino) != expected:
                raise ValueError("source changed during input preparation: " + str(path))
        backup = Path(tmp) / "previous"
        if out.exists():
            out.rename(backup)
        try:
            stage.rename(out)
        except BaseException:
            if backup.exists():
                backup.rename(out)
            raise
    print(f"PhenoRadar inputs: {out} ({len(species)} species)", flush=True)
    return dict(output=str(out), species=len(species), sections=inventory["sections"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True)
    parser.add_argument("--outdir")
    parser.add_argument("--exclude-species", default="[]", help="JSON exact species IDs")
    parser.add_argument("--trait", default="C4")
    args = vars(parser.parse_args())
    args["exclusions"] = json.loads(args.pop("exclude_species"))
    export(**args)
