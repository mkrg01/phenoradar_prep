#!/usr/bin/env python3
"""Export a manually selected species subset from completed results only."""
import argparse
from collections import Counter
import csv
from functools import lru_cache
import json
from pathlib import Path
import re
import shutil
import sqlite3
import tempfile

from common import file_record, now, write_json, write_tsv


SAFE = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]*\Z")
BUNDLES = {
    "odb": ["odb/merged/mappings.sqlite", "odb/merged/gene_orthogroups.tsv", "odb/merged/merge_qc.json"],
    "tpm": [f"tpm/{n}.tsv" for n in ["tpm", "tpm_wide", "tpm_sum", "tpm_sum_wide", "mapping_qc"]],
    "alignments": ["alignments/members.tsv", "alignments/provenance.json"],
    "kegg": [f"kegg/{n}.tsv" for n in ["genes", "gene_kos", "ko_tpm_sum", "ko_tpm_sum_wide", "ko_support", "mapping_qc"]],
    "phylogeny": ["phylogeny/species_tree.nwk", "phylogeny/species_tree.json", "phylogeny/gene_trees.nwk",
                  "phylogeny/gene_trees.json", "phylogeny/species_coverage.tsv"],
}
CONTRAST_BRANCHES = {"phylogeny": "metadata/samples.tsv",
                     "phylogeny_phenotyped": "phylogeny_phenotyped/selection/samples.tsv"}


def safe_name(value):
    if not isinstance(value, str) or not SAFE.fullmatch(value):
        raise ValueError(f"invalid identifier: {value!r}")
    return value


def validate_exclusions(values):
    if not isinstance(values, list):
        raise ValueError("exclude_species must be a YAML list of exact species IDs")
    for value in values:
        safe_name(value)
    if len(set(values)) != len(values):
        raise ValueError("exclude_species must not contain duplicate species IDs")
    return sorted(values)


def table(path, required=()):
    """Yield the header followed by validated rows, without loading large tables."""
    with Path(path).open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames
        if not fields or len(fields) != len(set(fields)) or not set(required) <= set(fields):
            raise ValueError(f"missing/duplicate columns: {path}; required {required}")
        yield fields
        for row in reader:
            if None in row or any(v is None for v in row.values()):
                raise ValueError(f"malformed TSV row: {path}:{reader.line_num}")
            yield row


def manifest(path):
    reader = table(path, ["species", "scientific_name", "run", "odb_species", "taxid"])
    fields = next(reader)
    rows, runs, species, filenames = [], {}, {}, {}
    for row in reader:
        for key in ["species", "run", "odb_species"]:
            safe_name(row[key])
        if row["run"] in runs:
            raise ValueError("duplicate run in sample manifest: " + row["run"])
        if row["species"] in species and any(row[k] != species[row["species"]][k]
                                              for k in ["scientific_name", "taxid", "odb_species"]):
            raise ValueError("conflicting metadata for species: " + row["species"])
        if row["odb_species"] in filenames and filenames[row["odb_species"]] != row["species"]:
            raise ValueError("colliding protein filenames")
        rows.append(row)
        runs[row["run"]] = row["species"]
        species[row["species"]] = row
        filenames[row["odb_species"]] = row["species"]
    if not rows:
        raise ValueError("empty source sample manifest")
    return fields, rows, species, runs


def discover(source, traits=None):
    """List external snapshot inputs; never request missing inference jobs."""
    source = Path(source).resolve()
    previous = source / "manifest.json"
    if previous.is_file() and json.loads(previous.read_text()).get("report_type") == "species_filter":
        raise ValueError("use the original analysis as source, not a previous filtered export")
    fields, rows, species, runs = manifest(source / "metadata/samples.tsv")
    files = {source / "metadata/samples.tsv"}
    sections = {}
    for name, names in BUNDLES.items():
        missing = [n for n in names if not (source / n).is_file()]
        sections[name] = {"status": "ready" if not missing else "absent" if len(missing) == len(names) else "incomplete",
                          "missing": missing}
        files.update(source / n for n in names if (source / n).is_file())
    for relative in ["run.json", "metadata/metadata_all.tsv", "metadata/metadata_high_busco.tsv", "metadata/selection.json",
                     "kegg/ko_modules.tsv", "kegg/ko_pathways.tsv", "kegg/reference_qc.json"]:
        if (source / relative).is_file():
            files.add(source / relative)
    proteins = [source / "proteins" / f"{r['odb_species']}_protein.fa" for r in species.values()]
    existing_proteins = [p for p in proteins if p.is_file()]
    sections["proteins"] = {"status": "ready" if len(existing_proteins) == len(proteins) else "incomplete" if existing_proteins else "absent",
                            "files_present": len(existing_proteins), "files_expected": len(proteins)}
    files.update(existing_proteins)
    if sections["alignments"]["status"] == "ready":
        report = json.loads((source / "alignments/provenance.json").read_text())
        for record in report["alignments"]:
            files.add(source / "alignments" / f"{safe_name(record['orthogroup'])}.faa")
    if sections["phylogeny"]["status"] == "ready":
        report = json.loads((source / "phylogeny/gene_trees.json").read_text())
        for record in report["retained"]:
            marker = safe_name(record["marker"])
            for folder in ["alignments", "alignments/raw"]:
                for suffix in [".faa", ".columns.tsv"]:
                    path = source / "phylogeny" / folder / (marker + suffix)
                    if path.is_file():
                        files.add(path)
        path = source / "phylogeny/dating/species_tree.dated.nwk"
        if path.is_file():
            files.add(path)
    if traits and Path(traits).is_file():
        files.add(Path(traits).resolve())
        sections["traits"] = {"status": "ready", "path": str(Path(traits).resolve())}
    else:
        sections["traits"] = {"status": "absent"}
    # Pair postprocessing needs only a completed rooted species tree, not gene
    # trees, prior pairs, or the inference programs. Discover both saved runs.
    metadata = source / "metadata/metadata_high_busco.tsv"
    score_columns = set(next(table(metadata))) if metadata.is_file() else set()
    for branch, samples in CONTRAST_BRANCHES.items():
        required = [f"{branch}/species_tree.nwk", f"{branch}/species_tree.json", samples,
                    "metadata/metadata_high_busco.tsv"]
        missing = [p for p in required if not (source / p).is_file()]
        if sections["traits"]["status"] != "ready":
            missing.append("species trait table")
        if not {"species", "busco_percent"} <= score_columns:
            missing.append("metadata columns species/busco_percent")
        present = any((source / p).is_file() for p in required[:2])
        sections[f"{branch}_contrast"] = {"status": "ready" if not missing else "incomplete" if present else "absent",
                                           "missing": missing}
        files.update(source / p for p in required if (source / p).is_file())
    return {"sections": sections, "files": sorted(str(p.resolve()) for p in files)}


def fasta_records(path):
    """Preserve full headers and residue case; wrap differences are cosmetic."""
    with Path(path).open() as handle:
        header, parts = None, []
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    if not parts:
                        raise ValueError(f"empty sequence: {path}")
                    yield header, "".join(parts)
                header, parts = line[1:], []
                if not header.split():
                    raise ValueError(f"empty FASTA identifier: {path}")
            elif header is None:
                raise ValueError(f"sequence before header: {path}")
            else:
                parts.append(line)
        if header is not None:
            if not parts:
                raise ValueError(f"empty sequence: {path}")
            yield header, "".join(parts)


class Export:
    def __init__(self, source, stage, exclusions, inventory):
        self.source, self.stage = source, stage
        self.fields, self.rows, by_species, self.runs = manifest(source / "metadata/samples.tsv")
        self.by_species, self.species = by_species, set(by_species)
        self.excluded = set(exclusions)
        if self.excluded - self.species:
            raise ValueError("exclude_species absent from selected source species: " + ", ".join(sorted(self.excluded - self.species)))
        self.keep = self.species - self.excluded
        if not self.keep:
            raise ValueError("exclude_species would remove every selected species")
        self.keep_runs = {r for r, s in self.runs.items() if s in self.keep}
        self.inventory, self.records, self.stats, self.counts = inventory, {}, {}, {}
        self.allowed_inputs = set(inventory["files"])
        self.ignored = ["phylogeny_phenotyped inference outputs (contrast pairs are recomputed)",
                        "contrast", "taxonomy_audit", "phylogeny/taxonomy_audit",
                        "raw logs, chunk results, and per-run computation caches"]

    def input(self, path, expected_hash=None):
        path = Path(path).resolve()
        key = str(path)
        if key not in self.allowed_inputs:
            raise ValueError(f"input not in frozen export inventory: {path}")
        if key not in self.records:
            stat = path.stat()
            self.stats[key] = (stat.st_size, stat.st_mtime_ns)
            self.records[key] = file_record(path)
        if expected_hash and self.records[key]["sha256"] != expected_hash:
            raise ValueError(f"source differs from its recorded checksum: {path}")
        return path

    def copy(self, source, destination):
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(self.input(source), destination)

    def run_row(self, row):
        if row["run"] not in self.runs or self.runs[row["run"]] != row["species"]:
            raise ValueError("table run/species differs from sample manifest: " + row["run"])
        return row["run"] in self.keep_runs

    def subset_table(self, relative, required, predicate, unique_runs=False, complete_runs=False):
        path = self.input(self.source / relative)
        reader = table(path, required)
        fields = next(reader)
        before = after = 0
        seen = set()
        dest = self.stage / relative
        dest.parent.mkdir(parents=True, exist_ok=True)
        with dest.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for row in reader:
                before += 1
                if "run" in row:
                    if unique_runs and row["run"] in seen:
                        raise ValueError(f"duplicate run in {relative}")
                    seen.add(row["run"])
                if predicate(row):
                    writer.writerow(row)
                    after += 1
        if complete_runs and seen != set(self.runs):
            raise ValueError(f"missing/extra runs in {relative}")
        self.counts[relative] = dict(before=before, after=after)

    def metadata(self):
        self.input(self.source / "metadata/samples.tsv")
        retained = [r for r in self.rows if r["species"] in self.keep]
        write_tsv(self.stage / "metadata/samples.tsv", self.fields, retained)
        write_tsv(self.stage / "excluded_samples.tsv", ["species", "run", "scientific_name", "taxid"],
                  [{k: r[k] for k in ["species", "run", "scientific_name", "taxid"]}
                   for r in self.rows if r["species"] in self.excluded])
        (self.stage / "metadata/species.txt").write_text("\n".join(sorted(self.keep)) + "\n")
        for filename in ["metadata_all.tsv", "metadata_high_busco.tsv"]:
            if (self.source / "metadata" / filename).is_file():
                # metadata_all also contains species removed by the original BUSCO
                # selection. The curated export contains only active sample rows.
                self.subset_table("metadata/" + filename, ["species", "run"],
                                  lambda row: row["run"] in self.keep_runs and self.run_row(row))
        if self.inventory["sections"]["traits"]["status"] == "ready":
            reader = table(self.input(self.inventory["sections"]["traits"]["path"]), ["species"])
            fields = next(reader)
            seen, rows = set(), []
            for row in reader:
                name = safe_name(row["species"].strip().replace(" ", "_"))
                if name in seen:
                    raise ValueError("duplicate normalized species in phenotype table: " + name)
                seen.add(name)
                if name in self.keep:
                    rows.append({**row, "species": name})
            write_tsv(self.stage / "metadata/species_trait.tsv", fields, rows)
        self.counts["species"] = dict(before=len(self.species), after=len(self.keep))
        self.counts["runs"] = dict(before=len(self.rows), after=len(retained))

    def proteins(self):
        destination = self.stage / "proteins"
        destination.mkdir()
        for name in sorted(self.keep):
            path = self.source / "proteins" / f"{self.by_species[name]['odb_species']}_protein.fa"
            (destination / path.name).symlink_to(self.input(path))
        self.counts["protein_files"] = dict(before=len(self.species), after=len(self.keep))

    def odb(self):
        destination = self.stage / "odb/merged"
        destination.mkdir(parents=True)
        path = self.input(self.source / "odb/merged/mappings.sqlite")
        with sqlite3.connect(destination / "mappings.sqlite", uri=True) as db:
            db.execute("ATTACH DATABASE ? AS original", (path.as_uri() + "?mode=ro",))
            unknown = db.execute("SELECT DISTINCT species FROM original.genes").fetchall()
            if {s for s, in unknown} != self.species:
                raise ValueError("ODB database species differ from the source manifest")
            db.executescript("""
                CREATE TABLE genes(query TEXT PRIMARY KEY, species TEXT NOT NULL);
                CREATE TABLE mappings(query TEXT NOT NULL REFERENCES genes(query), og TEXT NOT NULL,
                                      PRIMARY KEY(query, og)) WITHOUT ROWID;
                CREATE TEMP TABLE kept_species(species TEXT PRIMARY KEY);
                PRAGMA foreign_keys=ON;
            """)
            db.executemany("INSERT INTO kept_species VALUES (?)", ((s,) for s in sorted(self.keep)))
            db.execute("INSERT INTO genes SELECT g.query,g.species FROM original.genes g JOIN kept_species USING(species)")
            db.execute("CREATE INDEX genes_species ON genes(species)")
            db.execute("INSERT INTO mappings SELECT m.query,m.og FROM original.mappings m JOIN genes g USING(query)")
            if db.execute("SELECT 1 FROM original.mappings m LEFT JOIN original.genes g USING(query) WHERE g.query IS NULL LIMIT 1").fetchone():
                raise ValueError("orphan gene in source ODB mappings")
            pairs_before = db.execute("SELECT count(*) FROM original.mappings").fetchone()[0]
            pairs_after = db.execute("SELECT count(*) FROM mappings").fetchone()[0]
            genes_before = db.execute("SELECT count(*) FROM original.genes").fetchone()[0]
            genes_after = db.execute("SELECT count(*) FROM genes").fetchone()[0]
            write_tsv(destination / "gene_orthogroups.tsv", ["#query", "ODB_OG"],
                      ({"#query": g, "ODB_OG": og} for g, og in db.execute("SELECT query,og FROM mappings ORDER BY query,og")))
            # All copies and ambiguous assignments of retained genes remain intact.
            self.counts["odb_gene_og_pairs"] = dict(before=pairs_before, after=pairs_after)
            self.counts["odb_genes"] = dict(before=genes_before, after=genes_after)
            write_json(destination / "filter_qc.json", {"operation": "subset_original_mappings", "remapped": False,
                       "unique_gene_og_pairs": pairs_after, "genes": genes_after})

    def tpm(self):
        for relative in BUNDLES["tpm"]:
            unique = relative.endswith("_wide.tsv") or relative.endswith("mapping_qc.tsv")
            self.subset_table(relative, ["species", "run"], self.run_row, unique_runs=unique, complete_runs=True)

    def verify_gene_owners(self, db, table_name, gene_column):
        if self.inventory["sections"]["odb"]["status"] != "ready":
            return
        original = self.input(self.source / "odb/merged/mappings.sqlite")
        db.execute("ATTACH DATABASE ? AS original_odb", (original.as_uri() + "?mode=ro",))
        # Identifiers below are fixed internal table/column names, never config.
        bad = db.execute(f"SELECT m.{gene_column} FROM {table_name} m LEFT JOIN original_odb.genes g "
                         f"ON m.{gene_column}=g.query WHERE g.query IS NULL OR m.species!=g.species LIMIT 1").fetchone()
        if bad:
            raise ValueError("gene/species ownership differs from ODB: " + bad[0])

    def alignments(self):
        destination = self.stage / "alignments"
        destination.mkdir()
        report = json.loads(self.input(self.source / "alignments/provenance.json").read_text())
        member_path = self.input(self.source / "alignments/members.tsv", report.get("members", {}).get("sha256"))
        with sqlite3.connect(self.stage / ".alignment_members.sqlite", uri=True) as db:
            db.execute("CREATE TABLE members(og TEXT, gene TEXT, species TEXT, PRIMARY KEY(og,gene)) WITHOUT ROWID")
            reader = table(member_path, ["orthogroup", "gene_id", "species"])
            next(reader)
            batch = []
            for row in reader:
                if row["species"] not in self.species:
                    raise ValueError("unknown species in alignment membership")
                batch.append((safe_name(row["orthogroup"]), row["gene_id"], row["species"]))
                if len(batch) == 10000:
                    db.executemany("INSERT INTO members VALUES (?,?,?)", batch)
                    batch.clear()
            db.executemany("INSERT INTO members VALUES (?,?,?)", batch)
            self.verify_gene_owners(db, "members", "gene")
            if self.inventory["sections"]["odb"]["status"] == "ready":
                missing = db.execute("SELECT m.gene FROM members m LEFT JOIN original_odb.mappings p "
                                     "ON m.gene=p.query AND m.og=p.og WHERE p.query IS NULL LIMIT 1").fetchone()
                extra = db.execute("SELECT p.query FROM original_odb.mappings p LEFT JOIN members m "
                                   "ON m.gene=p.query AND m.og=p.og WHERE m.gene IS NULL LIMIT 1").fetchone()
                if missing or extra:
                    raise ValueError("alignment gene/OG membership differs from completed ODB mappings")
            declared = [safe_name(r["orthogroup"]) for r in report["alignments"]]
            if len(set(declared)) != len(declared) or set(declared) != {og for og, in db.execute("SELECT DISTINCT og FROM members")}:
                raise ValueError("alignment membership differs from completed alignment inventory")
            records, empty = [], []
            with (destination / "members.tsv").open("w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(["orthogroup", "gene_id", "species"])
                for record in report["alignments"]:
                    og = record["orthogroup"]
                    members = dict(db.execute("SELECT gene,species FROM members WHERE og=?", (og,)))
                    path = self.input(self.source / "alignments" / f"{og}.faa", record["alignment"]["sha256"])
                    seen, width, retained = set(), None, []
                    for header, sequence in fasta_records(path):
                        gene = header.split()[0]
                        if gene not in members or gene in seen:
                            raise ValueError(f"alignment gene differs from membership: {og}: {gene}")
                        seen.add(gene)
                        if width is None:
                            width = len(sequence)
                        if len(sequence) != width:
                            raise ValueError(f"unequal alignment lengths: {og}")
                        if members[gene] in self.keep:
                            retained.append((header, sequence))
                            writer.writerow([og, gene, members[gene]])
                    if seen != set(members):
                        raise ValueError("missing genes in alignment: " + og)
                    if retained:
                        with (destination / f"{og}.faa").open("w") as fasta:
                            for header, sequence in retained:
                                fasta.write(f">{header}\n{sequence}\n")
                    else:
                        empty.append(og)
                    records.append(dict(orthogroup=og, before=len(members), after=len(retained), columns=width))
            write_json(destination / "filter_qc.json", {"realigned": False, "columns_changed": False,
                       "empty_orthogroups": empty, "alignments": records})
            self.counts["alignment_sequences"] = dict(before=sum(r["before"] for r in records), after=sum(r["after"] for r in records))
        (self.stage / ".alignment_members.sqlite").unlink()

    def kegg(self):
        destination = self.stage / "kegg"
        destination.mkdir()
        with sqlite3.connect(self.stage / ".ko_genes.sqlite", uri=True) as db:
            db.execute("CREATE TABLE genes(gene TEXT PRIMARY KEY,species TEXT)")
            def gene_row(row):
                if row["species"] not in self.species:
                    raise ValueError("unknown species in KO genes")
                db.execute("INSERT INTO genes VALUES (?,?)", (row["gene_id"], row["species"]))
                return row["species"] in self.keep
            self.subset_table("kegg/genes.tsv", ["species", "gene_id"], gene_row)
            self.verify_gene_owners(db, "genes", "gene")
            @lru_cache(maxsize=10000)
            def owner(gene):
                value = db.execute("SELECT species FROM genes WHERE gene=?", (gene,)).fetchone()
                return value[0] if value else None
            def hit_row(row):
                if owner(row["gene_id"]) != row["species"]:
                    raise ValueError("KO gene/species differs from its membership table")
                return row["species"] in self.keep
            self.subset_table("kegg/gene_kos.tsv", ["species", "gene_id", "ko"], hit_row)
            for relative in BUNDLES["kegg"][2:]:
                unique = relative.endswith("_wide.tsv") or relative.endswith("mapping_qc.tsv")
                self.subset_table(relative, ["species", "run"], self.run_row, unique_runs=unique, complete_runs=unique)
            # These maps describe KO features, not species. Preserve source feature
            # axes, including zeros and missing values, across the table subset.
            for name in ["ko_modules.tsv", "ko_pathways.tsv", "reference_qc.json"]:
                path = self.source / "kegg" / name
                if path.is_file():
                    self.copy(path, destination / name)
        (self.stage / ".ko_genes.sqlite").unlink()


def export_contrast(job, branch, trait, seed):
    from contrast_pairs import from_tree
    from plot_contrast_tree import plot
    samples = job.input(job.source / CONTRAST_BRANCHES[branch])
    if {r["species"] for r in list(table(samples))[1:]} - job.species:
        raise ValueError("contrast species manifest contains species outside the source dataset")
    out = job.stage / branch / "contrast"
    report = from_tree(job.input(job.source / branch / "species_tree.nwk"),
                       job.input(job.source / branch / "species_tree.json"), samples,
                       job.input(job.source / "metadata/metadata_high_busco.tsv"),
                       job.input(job.inventory["sections"]["traits"]["path"]),
                       out, trait=trait, seed=seed, exclude_species=sorted(job.excluded))
    plot(out / "summary_tree.nwk", out / "species_metadata.tsv", out / "summary.json", out)
    job.counts[f"{branch}_contrast"] = {k: report[k] for k in
                                       ["source_species", "retained_species", "observed_species", "contrast_pairs"]}


def export(source, exclusions, outdir=None, traits=None, contrast_trait="C4", seed=12345):
    exclusions = validate_exclusions(exclusions)
    source = Path(source).resolve()
    requested_out = Path(outdir) if outdir else source / "filtered"
    if requested_out.is_symlink():
        raise ValueError("filtered output directory must not be a symlink")
    out = requested_out.resolve()
    inventory = discover(source, traits)
    if out == source or any(out == Path(p) or out in Path(p).parents for p in inventory["files"]):
        raise ValueError("filtered output directory overlaps original inputs")
    if out.exists() and any(out.iterdir()):
        marker = out / "manifest.json"
        if not marker.is_file() or json.loads(marker.read_text()).get("report_type") != "species_filter":
            raise ValueError("refusing to replace a directory not owned by species_filter")
    out.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".species-filter-", dir=out.parent) as tmp:
        stage = Path(tmp) / "export"
        stage.mkdir()
        job = Export(source, stage, exclusions, inventory)
        print(f"Species filter: {len(job.species)} -> {len(job.keep)} species", flush=True)
        job.metadata()
        for name, section in inventory["sections"].items():
            if name == "traits":
                continue
            print(f"{name}: {section['status']}", flush=True)
            if section["status"] == "ready":
                if name.endswith("_contrast"):
                    export_contrast(job, name.removesuffix("_contrast"), contrast_trait, seed)
                elif name == "phylogeny":
                    from filter_species_phylogeny import export_phylogeny
                    export_phylogeny(job)
                else:
                    getattr(job, name)()
        # Record all known snapshot metadata, including the sources used to
        # decide that an optional branch is incomplete. Never edit those files.
        for path in inventory["files"]:
            job.input(path)
        for path, expected in job.stats.items():
            stat = Path(path).stat()
            if (stat.st_size, stat.st_mtime_ns) != expected:
                raise ValueError("source changed during export: " + path)
        outputs = []
        for path in sorted(stage.rglob("*")):
            if path.is_file():
                if path.is_symlink():
                    record = dict(job.records[str(path.resolve())], storage="symlink", target=str(path.resolve()))
                else:
                    record = dict(file_record(path), storage="file")
                outputs.append({**record, "path": str(path.relative_to(stage))})
        summary = {"report_type": "species_filter", "schema_version": 2, "created_at": now(),
                   "source": str(source), "exclude_species": exclusions, "retained_species": sorted(job.keep),
                   "retained_runs": sorted(job.keep_runs), "counts": job.counts, "sections": inventory["sections"],
                   "inputs": list(job.records.values()), "outputs": outputs, "not_exported": job.ignored,
                   "contrast": {"trait": contrast_trait, "seed": seed},
                   "code": [file_record(Path(__file__).with_name(name)) for name in
                            ["filter_species.py", "filter_species_phylogeny.py", "contrast_pairs.py",
                             "plot_contrast_tree.py", "species_traits.py", "phylogeny_root.py", "common.py"]],
                   "policies": {"sequence_inference_recomputed": False,
                                "contrast_pairs": "recomputed from completed molecular trees after exclusions; NCBI representative analysis untouched",
                                "expression_values": "original strings preserved",
                                "feature_columns": "source OG/KO axes retained, including all-zero and unavailable columns",
                                "alignment_columns": "unchanged, including columns left all-gap",
                                "protein_files": "symlinks to retained species files; keep original results available",
                                "phylogenies": "pruned derivatives only; no inference/support/calibration refitting"}}
        write_json(stage / "manifest.json", summary)
        # Restore the previous complete export if publication itself fails.
        backup = Path(tmp) / "previous"
        if out.exists():
            out.rename(backup)
        try:
            stage.rename(out)
        except BaseException:
            if backup.exists():
                backup.rename(out)
            raise
    print(json.dumps({"output": str(out), "species": len(job.keep), "runs": len(job.keep_runs)}), flush=True)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True)
    parser.add_argument("--outdir")
    parser.add_argument("--traits")
    parser.add_argument("--contrast-trait", default="C4")
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument("--exclude-species", default="[]", help="JSON list of exact species IDs")
    args = vars(parser.parse_args())
    args["exclusions"] = json.loads(args.pop("exclude_species"))
    export(**args)


if __name__ == "__main__":
    main()
