#!/usr/bin/env python3
"""Publish only PhenoRadar inputs from completed analysis outputs."""
import argparse
import csv
import gzip
import json
import math
from pathlib import Path
import tempfile

from common import sha256, species_from_gene_id, write_tsv
from filter_species import fasta_records, manifest, safe_name, table, validate_exclusions
from phenoradar_metadata import read_base, with_pairs
from layout import ORTHOGROUP_EXPRESSION, ORTHOGROUP_ALIGNMENTS, CONTRAST_BRANCHES


DEFAULTS = dict(orthogroups="auto", kegg="auto", alignments="auto",
                kegg_groups=["module"], orthogroup_annotations=None, contrast=None, tree=None)
CONTRASTS = CONTRAST_BRANCHES


def validate_settings(settings=None):
    if settings is not None and not isinstance(settings, dict):
        raise ValueError("phenoradar must be a configuration mapping")
    settings = settings or {}
    if set(settings) - set(DEFAULTS):
        raise ValueError(f"unknown phenoradar settings: {sorted(set(settings) - set(DEFAULTS))}")
    result = {**DEFAULTS, **settings}
    for key in ["orthogroups", "kegg", "alignments"]:
        if type(result[key]) is not bool and result[key] != "auto":
            raise ValueError(f"phenoradar.{key} must be true, false, or auto")
    groups = result["kegg_groups"]
    if (not isinstance(groups, list) or any(not isinstance(g, str) or g not in {"module", "pathway"} for g in groups)
            or len(set(groups)) != len(groups)):
        raise ValueError("phenoradar.kegg_groups must list module and/or pathway, without duplicates")
    if result["contrast"] is not None and (not isinstance(result["contrast"], str) or result["contrast"] not in CONTRASTS):
        raise ValueError("phenoradar.contrast must be null or one of " + ", ".join(sorted(CONTRASTS)))
    for key in ["orthogroup_annotations", "tree"]:
        value = result[key]
        if value is not None and (not isinstance(value, str) or not value.strip()):
            raise ValueError(f"phenoradar.{key} must be null or a file path")
    return result


def discover(source, settings=None, exclusions=()):
    """Freeze existing files as absolute inputs, without requesting producers."""
    settings = validate_settings(settings)
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
            raise ValueError(f"missing {path}; run phenoradar_metadata (then filter_species if excluding species)")
        files.add(path)

    sections, links, alignment_hashes = {}, {}, {}

    def branch(name, required, published):
        mode = settings[name]
        if mode is False:
            sections[name] = "disabled"
            return False
        missing = [p for p in required if not (source / p).is_file()]
        present = any((source / p).exists() for p in required)
        status = "ready" if not missing else "incomplete" if present else "absent"
        sections[name] = status
        files.update(source / p for p in required if (source / p).is_file())
        if missing:
            if mode is True:
                raise ValueError(f"required {name} inputs are {status}: {', '.join(missing)}")
            return False
        links.update({dest: source / relative for dest, relative in published.items()})
        return True

    branch("orthogroups", [f"{ORTHOGROUP_EXPRESSION}/tpm.tsv", f"{ORTHOGROUP_EXPRESSION}/mapping_qc.tsv"], {"tpm.tsv": f"{ORTHOGROUP_EXPRESSION}/tpm.tsv"})
    maps = {f"kegg/ko_{group}s.tsv": f"kegg/ko_{group}s.tsv" for group in settings["kegg_groups"]}
    branch("kegg", ["kegg/ko_tpm_sum.tsv", "kegg/mapping_qc.tsv", "kegg/ko_support.tsv",
                    "kegg/reference_qc.json", *maps],
           {"kegg/ko_tpm_sum.tsv": "kegg/ko_tpm_sum.tsv", **maps})
    completion = f"{ORTHOGROUP_ALIGNMENTS}/filter_qc.json" if filtered else f"{ORTHOGROUP_ALIGNMENTS}/provenance.json"
    if branch("alignments", [completion], {}):
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
            links[f"alignments/{og}.faa"] = path
            files.add(path)
            if not filtered:
                alignment_hashes[str(path)] = record["alignment"]["sha256"]
    elif sections["alignments"] == "absent" and any((source / ORTHOGROUP_ALIGNMENTS).glob("*.faa")):
        sections["alignments"] = "incomplete"

    if not ({"tpm.tsv", "kegg/ko_tpm_sum.tsv"} & links.keys()):
        raise ValueError("no completed expression input selected; run the OG TPM or kegg step first")
    pairs = None
    if settings["contrast"]:
        pairs = source / settings["contrast"] / "species_metadata.tsv"
        for path in [pairs, pairs.with_name("summary.json")]:
            if not path.is_file():
                raise ValueError("selected contrast output is incomplete: " + str(path))
            files.add(path)
    for key, destination in [("tree", "species_tree.nwk"),
                             ("orthogroup_annotations", "orthogroup_annotations.tsv")]:
        if settings[key]:
            path = Path(settings[key]).resolve()
            if not path.is_file():
                raise ValueError(f"selected {key} is missing: {path}")
            if key == "orthogroup_annotations" and path.suffix == ".gz":
                destination += ".gz"
            links[destination] = path
            files.add(path)
    return dict(source=str(source), files=sorted(str(p) for p in files), links=links,
                sections=sections, filtered=filtered, pairs=pairs, alignment_hashes=alignment_hashes,
                original_runs={r["run"]: r["species"] for r in original_samples if r["species"] not in exclusions})


def validate_expression(path, species_runs, feature, value):
    """Check producer long tables using memory bounded by one species' features."""
    reader = table(path, ["species", "run", feature, value])
    next(reader)
    species_seen, features, current_features = set(), set(), set()
    current = None
    for row in reader:
        species = row["species"]
        if species_runs.get(species) != row["run"]:
            raise ValueError(f"expression run/species differs from samples: {path}: {species}/{row['run']}")
        if species != current:
            if species in species_seen:
                raise ValueError(f"duplicate/discontiguous species block in expression: {path}: {species}")
            species_seen.add(species)
            current, current_features = species, set()
        name = safe_name(row[feature])
        if name in current_features:
            raise ValueError(f"duplicate species/feature coordinate: {path}: {species}/{name}")
        current_features.add(name)
        features.add(name)
        try:
            number = float(row[value])
        except ValueError as error:
            raise ValueError(f"invalid numeric expression: {path}: {row[value]!r}") from error
        if not math.isfinite(number) or number < 0:
            raise ValueError(f"expression must be finite and nonnegative: {path}")
    if species_seen != set(species_runs):
        raise ValueError(f"expression lacks selected species: {path}: {sorted(set(species_runs) - species_seen)}")
    return features


def validate_qc(path, species_runs):
    reader = table(path, ["species", "run"])
    next(reader)
    seen = set()
    for row in reader:
        species = row["species"]
        if species in seen or species_runs.get(species) != row["run"]:
            raise ValueError("expression QC run/species differs from samples: " + str(path))
        seen.add(species)
    if seen != set(species_runs):
        raise ValueError("expression QC lacks selected species: " + str(path))


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


def validate_tree(path, species):
    from ete4 import Tree
    text = path.read_text().strip()
    if text.count(";") != 1:
        raise ValueError("expected one Newick species tree")
    tree = Tree(text, parser=1)
    names = list(tree.leaf_names())
    if len(names) != len(set(names)) or set(names) != species:
        raise ValueError("tree species differ from selected metadata; use the matching full/pruned tree")
    for node in tree.traverse():
        if node.dist is not None and (not math.isfinite(node.dist) or node.dist < 0):
            raise ValueError("tree branch lengths must be finite and nonnegative when present")


def validate_annotations(path, orthogroups):
    if not orthogroups:
        raise ValueError("OG annotations require an OG expression input")
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
    if not matched:
        raise ValueError("OG annotations do not match expression IDs; check the OrthoDB reference version")
    print(f"OG descriptions: {len(matched)}/{len(orthogroups)} expressed OGs", flush=True)


def check_destination(out):
    if out.is_symlink():
        raise ValueError("phenoradar_inputs output directory must not be a symlink")
    allowed = {"species_metadata.tsv", "tpm.tsv", "species_tree.nwk", "orthogroup_annotations.tsv",
               "orthogroup_annotations.tsv.gz", ".snakemake_timestamp"}
    if out.exists():
        for path in out.iterdir():
            if path.name in allowed and (path.is_file() or path.is_symlink()):
                continue
            if path.name in {"alignments", "kegg"} and path.is_dir() and not path.is_symlink():
                for child in path.iterdir():
                    valid = (child.suffix == ".faa" if path.name == "alignments" else
                             child.name in {"ko_tpm_sum.tsv", "ko_modules.tsv", "ko_pathways.tsv"})
                    if not valid or not child.is_symlink():
                        raise ValueError("unrecognized file in output directory: " + str(child))
                continue
            raise ValueError("refusing to replace unrecognized output contents: " + str(path))


def export(source, outdir=None, settings=None, exclusions=(), trait="C4"):
    out = Path(outdir) if outdir else Path(source) / "phenoradar_inputs"
    check_destination(out)
    out = out.resolve()
    inventory = discover(source, settings, exclusions)
    source = Path(inventory["source"])
    paths = [Path(p) for p in inventory["files"]]
    if out == source or any(out == p.resolve() or out in p.resolve().parents for p in paths):
        raise ValueError("phenoradar_inputs output overlaps its source files")
    stats = {p: (p.stat().st_size, p.stat().st_mtime_ns, p.stat().st_ino) for p in paths}
    print("PhenoRadar source: " + str(source), flush=True)
    for name, status in inventory["sections"].items():
        print(f"{name}: {status}", flush=True)
    expected_hashes = dict(inventory["alignment_hashes"])
    if inventory["sections"]["kegg"] == "ready":
        reference = json.loads((source / "kegg/reference_qc.json").read_text())
        results = reference.get("results", [])
        maps = {Path(r["path"]).name: r["sha256"] for r in results}
        for dest, path in inventory["links"].items():
            if dest in {"kegg/ko_modules.tsv", "kegg/ko_pathways.tsv"}:
                if path.name not in maps:
                    raise ValueError("KO membership is not recorded in KEGG reference completion: " + str(path))
                expected_hashes[str(path)] = maps[path.name]
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
    _, samples, by_species, _ = manifest(source / "metadata/samples.tsv")
    species_runs = {}
    for row in samples:
        if row["species"] in species_runs:
            raise ValueError("PhenoRadar requires one run per species; explicitly select runs before export: " + row["species"])
        species_runs[row["species"]] = row["run"]
    species = set(by_species)
    if filtered and species != set(filtered["retained_species"]):
        raise ValueError("filtered samples differ from recorded retained species")
    if filtered and {r["run"]: r["species"] for r in samples} != inventory["original_runs"]:
        raise ValueError("filtered runs differ from the current original samples; rerun filter_species")
    metadata = read_base(source / "metadata/species_metadata.tsv", trait)
    if {r["species"] for r in metadata} != species:
        raise ValueError("PhenoRadar metadata species differ from samples")
    if inventory["pairs"]:
        metadata = with_pairs(metadata, inventory["pairs"], trait)
    links, orthogroups = inventory["links"], set()
    for dest, feature, value, qc in [("tpm.tsv", "orthogroup", "tpm", f"{ORTHOGROUP_EXPRESSION}/mapping_qc.tsv"),
                                    ("kegg/ko_tpm_sum.tsv", "ko", "tpm_sum", "kegg/mapping_qc.tsv")]:
        if dest in links:
            print("Checking expression: " + str(links[dest]), flush=True)
            features = validate_expression(links[dest], species_runs, feature, value)
            validate_qc(source / qc, species_runs)
            if feature == "orthogroup":
                orthogroups = features
    for dest, path in links.items():
        if dest.startswith("alignments/"):
            validate_alignment(path, species)
        elif dest in {"kegg/ko_modules.tsv", "kegg/ko_pathways.tsv"}:
            validate_groups(path, "module" if "modules" in dest else "pathway")
        elif dest == "species_tree.nwk":
            validate_tree(path, species)
        elif dest.startswith("orthogroup_annotations."):
            validate_annotations(path, orthogroups)
    out.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".phenoradar-inputs-", dir=out.parent) as tmp:
        stage = Path(tmp) / "inputs"
        stage.mkdir()
        if inventory["pairs"]:
            write_tsv(stage / "species_metadata.tsv", ["species", trait, "contrast_pair_id", "family"], metadata)
        else:
            (stage / "species_metadata.tsv").symlink_to(source / "metadata/species_metadata.tsv")
        for dest, path in links.items():
            target = stage / dest
            target.parent.mkdir(parents=True, exist_ok=True)
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
    parser.add_argument("--settings", default="{}", help="JSON phenoradar settings")
    parser.add_argument("--exclude-species", default="[]", help="JSON exact species IDs")
    parser.add_argument("--trait", default="C4")
    args = vars(parser.parse_args())
    args["settings"] = json.loads(args["settings"])
    args["exclusions"] = json.loads(args.pop("exclude_species"))
    export(**args)
