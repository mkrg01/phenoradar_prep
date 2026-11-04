#!/usr/bin/env bash
set -euo pipefail

# Self-contained so both Snakemake and the image builder hash this repair.
"${CONDA_PREFIX:?}/bin/python" - <<'PY'
import hashlib
import json
import os
from pathlib import Path
import re

prefix = Path(os.environ["CONDA_PREFIX"]).absolute()
python = prefix / "bin/python"
if not python.is_file() or not os.access(python, os.X_OK):
    raise RuntimeError(f"ODB environment Python is not executable: {python}")
if re.search(r"\s", str(python)):
    raise RuntimeError("Orthologer requires an environment path without whitespace")
root = prefix / "orthologer/bin"
if not (root / "check_fasta.py").is_file():
    raise RuntimeError(f"Orthologer FASTA checker is missing: {root}")
manifest = prefix / "share/phenoradar/odb-runtime.json"
previous = json.loads(manifest.read_text()) if manifest.exists() else {}
previous = {entry["path"]: entry for entry in previous.get("scripts", [])}
records = []
for path in sorted(root.glob("*.py")):
    content = path.read_bytes()
    first, separator, body = content.partition(b"\n")
    if not first.startswith(b"#!") or b"python" not in first:
        continue
    if path.is_symlink():
        raise RuntimeError(f"Refusing to rewrite a symlink: {path}")
    match = re.fullmatch(rb"#!\S*/(?:env\s+)?python(?:3(?:\.\d+)?)?(\s+.*)?", first)
    if not match:
        raise RuntimeError(f"Unsupported Orthologer interpreter: {first!r}")
    repaired = b"#!" + os.fsencode(python) + (match[1] or b"") + separator + body
    relative = str(path.relative_to(prefix))
    before = hashlib.sha256(content).hexdigest()
    after = hashlib.sha256(repaired).hexdigest()
    prior = previous.get(relative, {})
    record = {"path": relative, "original_sha256": before,
              "original_shebang": first.decode(), "sha256": after}
    if prior.get("sha256") == before:
        record.update({key: prior[key] for key in ("original_sha256", "original_shebang")})
    if repaired != content:
        path.write_bytes(repaired)  # Preserve the installed executable mode.
    records.append(record)
manifest.parent.mkdir(parents=True, exist_ok=True)
manifest.write_text(json.dumps({"python": str(python), "scripts": records}, indent=2) + "\n")
print(f"Verified {len(records)} Orthologer Python entry points: {python}")
PY
