"""ODB subprocess and reference helpers shared by setup and mapping."""
import os
import shutil
import subprocess
from pathlib import Path

from common import file_record


def odb_environment(prefix, version):
    env = os.environ.copy()
    # ODB exports API/work variables. Do not inherit another project's settings.
    for key in list(env):
        if key.startswith("ODBAPI_") or key in {"ODBMAPPER_WORK", "MAP_ORTHODB_DATA", "DBI_DOWNLOAD"}:
            env.pop(key)
    env["ODBAPI_URL_VERSION"] = version
    if prefix:
        prefix = str(Path(prefix).resolve())
        if not Path(prefix, "bin").is_dir():
            raise ValueError(f"ODB environment not found: {prefix}")
        env["PATH"] = str(Path(prefix, "bin")) + os.pathsep + env["PATH"]
        env["CONDA_PREFIX"] = prefix
    return env


def command_path(command, env):
    path = shutil.which(command, path=env["PATH"])
    if path is None:
        raise ValueError(f"ODB executable not found: {command}")
    return str(Path(path).absolute())


def software_records(command, prefix):
    records = [file_record(command)]
    if not prefix and (Path(command).parent.parent / "conda-meta").is_dir():
        prefix = str(Path(command).parent.parent)
    if prefix:
        # Package records pin every installed dependency; hash actual Orthologer
        # scripts too, including local modifications to the installed package.
        paths = list(Path(prefix, "conda-meta").glob("*.json"))
        root = Path(prefix, "orthologer")
        for directory in ["bin", "etc"]:
            paths.extend(p for p in (root / directory).rglob("*") if p.is_file())
        records.extend(file_record(p) for p in sorted(set(paths)))
    return records


def check_storage(path, min_free_gb, allow_nonlocal=False):
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    if min_free_gb < 0:
        raise ValueError("minimum free space cannot be negative")
    kind = subprocess.check_output(["findmnt", "-n", "-o", "FSTYPE", "-T", str(path)], text=True).strip()
    if not allow_nonlocal and kind not in {"ext2", "ext3", "ext4", "xfs", "btrfs"}:
        raise ValueError(f"ODB work is on {kind}; configure local storage or explicitly allow nonlocal filesystems")
    free = shutil.disk_usage(path).free
    if free < min_free_gb * 1024**3:
        raise ValueError(f"ODB needs at least {min_free_gb} GiB free at {path}; available {free / 1024**3:.1f} GiB")
