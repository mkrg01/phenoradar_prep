#!/usr/bin/env python3
"""Check native Orthologer FASTA import before scheduling mapping chunks."""
import argparse
import os
from pathlib import Path
import re
import subprocess
import tempfile

from common import file_record, now, write_json
from odb_environment import command_path, odb_environment

PROTEIN = ">query\nMKWVTFISLLFLFSSAYSRGVFRRDTHKSEIAHRFKDLGE\n"


def check(command="ODB-mapper", prefix="", version="v12", work_dir=None):
    env = odb_environment(prefix, version)
    command = Path(command_path(command, env))
    candidates = [Path(prefix)] if prefix else [command.parent.parent, command.resolve().parent.parent]
    installation = next((p for p in candidates if (p / "orthologer/bin/check_fasta.py").is_file()), None)
    if installation is None:
        raise ValueError(f"Orthologer FASTA checker not found next to {command}; specify its --prefix")
    installation = installation.absolute()
    helper = installation / "orthologer/bin/check_fasta.py"
    orthologer = installation / "bin/orthologer"
    env["PATH"] = str(installation / "bin") + os.pathsep + env.get("PATH", "")
    env["CONDA_PREFIX"] = str(installation)
    # Keep user configuration and temporary outputs away from real projects.
    for key in list(env):
        if key.startswith(("DIR_ORTHOLOGER", "DIR_PROJECT", "ORTHOLOGER_")):
            env.pop(key)
    with tempfile.TemporaryDirectory(prefix="odb-runtime-", dir=work_dir) as temporary:
        work = Path(temporary).resolve()
        if not re.fullmatch(r"[A-Za-z0-9_./+-]+", str(work)):
            raise ValueError(f"ODB path contains unsupported characters: {work}")
        env.update(HOME=str(work), TMPDIR=str(work), PYTHONDONTWRITEBYTECODE="1")
        fasta = work / "input.faa"
        fasta.write_text(PROTEIN)
        commands = []

        def execute(argv):
            try:
                result = subprocess.run(list(map(str, argv)), cwd=work, env=env, text=True,
                                        capture_output=True, timeout=60, check=True)
            except (OSError, subprocess.TimeoutExpired, subprocess.CalledProcessError) as error:
                outputs = [getattr(error, name, "") for name in ("stdout", "stderr")]
                detail = "\n".join(value.decode(errors="replace") if isinstance(value, bytes)
                                   else value or "" for value in outputs)
                raise RuntimeError(f"Orthologer runtime check failed: {argv[0]}\n{error}\n{detail[-8000:]}") from error
            commands.append({"argv": list(map(str, argv)), "returncode": result.returncode})
            return result.stdout.strip()

        # Execute directly: using `python check_fasta.py` would hide a broken shebang.
        detected = execute([helper, "-f", fasta, "--type"])
        if detected != "PROT":
            raise ValueError(f"Orthologer misclassified the smoke protein: {detected!r}")
        project = work / "project"
        project.mkdir()
        execute([orthologer, "-c", "create", "-w", project])
        manifest = work / "input.txt"
        manifest.write_text(f"+smoke {fasta}\n")
        execute([orthologer, "-c", "import", "-w", project, "-d", manifest, "-p", "smoke"])
        # Upstream may report exit 0 even when individual imports fail.
        sequences = list((project / "Sequences").glob("smoke.*"))
        if len(sequences) != 1 or sequences[0].read_text() != PROTEIN.replace(">query", ">smoke:000000"):
            raise ValueError("Orthologer did not import the expected protein sequence")
        return {"completed_at": now(), "checks": {"fasta_type": detected, "fasta_import": "passed"},
                "software": [file_record(command), file_record(helper), file_record(orthologer)],
                "commands": commands}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--command", default="ODB-mapper")
    parser.add_argument("--prefix", default="")
    parser.add_argument("--version", default="v12")
    parser.add_argument("--output", required=True)
    args = vars(parser.parse_args())
    output = Path(args.pop("output"))
    output.parent.mkdir(parents=True, exist_ok=True)
    write_json(output, check(**args, work_dir=output.parent))
