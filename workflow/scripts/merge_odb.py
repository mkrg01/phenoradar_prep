#!/usr/bin/env python3
"""Validate ODB mappings and build an indexed, deduplicated gene-to-OG table."""
import argparse
import json
import os
import sqlite3
import tempfile
from pathlib import Path

from common import file_record, now, read_tsv, write_json, write_tsv
from translate_cds import fasta_ids


def annotation_pairs(path):
    """ODB v12 emits either annotated hits or bare cluster hits (no descriptions)."""
    indices = None
    saw_empty = False
    with open(path) as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            fields = line.rstrip("\r\n").split("\t")
            if line.startswith("#"):
                if "#query" in fields and "ODB_OG" in fields:
                    indices = fields.index("#query"), fields.index("ODB_OG")
                elif "#cluster_id" in fields and "gene_id" in fields:
                    indices = fields.index("gene_id"), fields.index("#cluster_id")
                elif line.startswith("# No data found"):
                    saw_empty = True
                continue
            if indices is None or len(fields) <= max(indices):
                raise ValueError(f"unrecognized ODB table at {path}:{number}")
            query, og = (fields[i] for i in indices)
            if not query or not og or og in {"-", "NA", "null"}:
                raise ValueError(f"empty query/OG at {path}:{number}")
            yield query, og
    if indices is None and not saw_empty:
        raise ValueError(f"ODB table has no recognized header: {path}")


def merge(samples, chunks, chunk_dir, protein_dir, database, mappings, qc,
          existing=None, version="v12", node=3193, source_plan=None):
    species = {row["species"]: row["odb_species"] for row in read_tsv(samples)}
    if source_plan:
        plan = json.loads(Path(source_plan).read_text())
        if plan["version"] != version or plan["node"] != node or set(plan["proteins"]) != set(species):
            raise ValueError("ODB plan differs from selected inputs")
        members = plan["sources"]
        for name, odb_name in species.items():
            if file_record(Path(protein_dir) / f"{odb_name}_protein.fa")["sha256"] != plan["proteins"][name]["sha256"]:
                raise ValueError(f"protein changed after ODB planning: {name}")
    elif existing:
        members = [{"kind": "existing", "root": str(existing), "species": sorted(species)}]
    else:
        members = [{"kind": "mapped", **c} for c in json.loads(Path(chunks).read_text())]
    if sorted(s for c in members for s in c["species"]) != sorted(species):
        raise ValueError("chunk membership differs from selected species")
    target = Path(database)
    target.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(dir=target.parent, prefix=".mappings.")
    os.close(fd)
    db = sqlite3.connect(temporary)
    try:
        db.executescript("""
            CREATE TABLE genes (query TEXT PRIMARY KEY, species TEXT NOT NULL);
            CREATE TABLE mappings (query TEXT NOT NULL REFERENCES genes(query), og TEXT NOT NULL,
                                   PRIMARY KEY (query, og)) WITHOUT ROWID;
            CREATE TEMP TABLE chunk_mappings (query TEXT, og TEXT);
            CREATE TEMP TABLE allowed_species (species TEXT PRIMARY KEY);
            PRAGMA foreign_keys = ON;
        """)
        with db:
            for name, odb_name in sorted(species.items()):
                path = Path(protein_dir) / f"{odb_name}_protein.fa"
                try:
                    db.executemany("INSERT INTO genes VALUES (?, ?)", ((gene, name) for gene in fasta_ids(path)))
                except sqlite3.IntegrityError as error:
                    raise ValueError("FASTA query IDs are not globally unique across species") from error
            db.execute("CREATE INDEX genes_species ON genes(species)")
        input_rows = excluded_rows = 0
        sources = []
        for chunk in members:
            imported = chunk["kind"] == "existing"
            label = "existing" if imported else chunk["chunk"]
            root = Path(chunk["root"]) if imported else Path(chunk_dir) / label
            provenance = root / ("snapshot.json" if imported else "provenance.json")
            subset = False
            if imported:
                from incremental_odb import load_snapshot
                if chunk.get("snapshot_sha256") and file_record(provenance)["sha256"] != chunk["snapshot_sha256"]:
                    raise ValueError("ODB snapshot changed after planning")
                snapshot, original = load_snapshot(root, version, node)
                if set(chunk["species"]) - original.keys():
                    raise ValueError("existing ODB snapshot has duplicate or missing selected species")
                for name in chunk["species"]:
                    odb_name = species[name]
                    path = Path(protein_dir) / f"{odb_name}_protein.fa"
                    if original[name]["odb_species"] != odb_name or file_record(path)["sha256"] != original[name]["sha256"]:
                        raise ValueError(f"protein differs from existing ODB input: {name}; remapping is required")
                subset = set(chunk["species"]) != set(original)
                record = {"results": [snapshot["annotations"]]}
            else:
                record = json.loads(provenance.read_text())
            for entry in record["results"]:
                path = root / Path(entry["path"]).name
                if file_record(path)["sha256"] != entry["sha256"]:
                    raise ValueError(f"ODB result changed after completion: {path}")
            annotated = root / ("annotations.tsv" if imported else f"{label}.og.annotations")
            batch = []
            with db:
                db.execute("DELETE FROM chunk_mappings")
                db.execute("DELETE FROM allowed_species")
                db.executemany("INSERT INTO allowed_species VALUES (?)", ((s,) for s in chunk["species"]))
                for query, og in annotation_pairs(annotated):
                    batch.append((query, og))
                    input_rows += 1
                    if len(batch) >= 10000:
                        db.executemany("INSERT INTO chunk_mappings VALUES (?, ?)", batch)
                        batch.clear()
                db.executemany("INSERT INTO chunk_mappings VALUES (?, ?)", batch)
                invalid = "FROM chunk_mappings c LEFT JOIN genes g ON c.query=g.query LEFT JOIN allowed_species a ON g.species=a.species WHERE a.species IS NULL"
                unknown = db.execute("SELECT c.query " + invalid + " LIMIT 1").fetchone()
                if unknown and not subset:
                    raise ValueError(f"{label}: query does not belong to this chunk's FASTA inputs: {unknown[0]}")
                if subset:
                    excluded_rows += db.execute("SELECT count(*) " + invalid).fetchone()[0]
                db.execute("INSERT OR IGNORE INTO mappings SELECT c.query,c.og FROM chunk_mappings c "
                           "JOIN genes g ON c.query=g.query JOIN allowed_species a ON g.species=a.species")
            sources.append(file_record(provenance))
        count = db.execute("SELECT count(*) FROM mappings").fetchone()[0]
        ambiguous = db.execute("SELECT count(*) FROM (SELECT query FROM mappings GROUP BY query HAVING count(*) > 1)").fetchone()[0]
        write_tsv(mappings, ["#query", "ODB_OG"],
                  ({"#query": q, "ODB_OG": og} for q, og in db.execute("SELECT query, og FROM mappings ORDER BY query, og")))
        kinds = {c["kind"] for c in members}
        mode = "mixed" if len(kinds) > 1 else "existing" if kinds == {"existing"} else "mapped"
        write_json(qc, {"created_at": now(), "chunks": len(members), "input_annotation_rows": input_rows,
                        "unique_gene_og_pairs": count, "duplicate_pairs_removed": input_rows - excluded_rows - count,
                        "excluded_annotation_rows": excluded_rows, "mode": mode,
                        "reused_species": sum(len(c["species"]) for c in members if c["kind"] == "existing"),
                        "mapped_species": sum(len(c["species"]) for c in members if c["kind"] == "mapped"),
                        "ambiguous_genes": ambiguous, "sources": sources})
        db.close()
        os.replace(temporary, target)
    finally:
        db.close()
        if os.path.exists(temporary):
            os.unlink(temporary)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ["samples", "chunks", "chunk-dir", "protein-dir", "database", "mappings", "qc"]:
        parser.add_argument(f"--{flag}", required=True)
    parser.add_argument("--existing")
    parser.add_argument("--source-plan")
    parser.add_argument("--version", default="v12")
    parser.add_argument("--node", type=int, default=3193)
    merge(**vars(parser.parse_args()))
