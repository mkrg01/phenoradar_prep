"""Validate the shared launcher for direct execution and one Slurm allocation."""
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def batch_workspace(tmp_path):
    checkout = tmp_path / "checkout"
    (checkout / "workflow").mkdir(parents=True)
    (checkout / "workflow/Snakefile").touch()
    shutil.copy2(ROOT / "run_pipeline.sh", checkout / "run_pipeline.sh")
    spool = tmp_path / "spool"
    spool.mkdir()
    script = spool / "slurm_script"
    shutil.copy2(ROOT / "run_pipeline.sh", script)
    env = {k: v for k, v in os.environ.items()
           if not k.startswith(("SLURM_", "SBATCH_", "SRUN_"))}
    env.update({"SLURM_JOB_ID": "123", "SLURM_CPUS_PER_TASK": "2",
                "SLURM_MEM_PER_NODE": "8192", "SLURM_NTASKS": "1",
                "XDG_CACHE_HOME": str(tmp_path / "inherited_cache")})
    return checkout, script, env


@pytest.mark.parametrize("mode,exit_code", [("batch_node", 0), ("batch_cpu", 7),
                                         ("direct", 0), ("direct", 7)])
def test_launcher_arguments_and_exit_status(batch_workspace, tmp_path, mode, exit_code):
    checkout, script, env = batch_workspace
    capture = tmp_path / "capture.json"
    fake = tmp_path / "fake snakemake"
    fake.write_text(
        f"#!{sys.executable}\n"
        "import json, os, sys\n"
        "from pathlib import Path\n"
        "Path(os.environ['PHENORADAR_TEST_CAPTURE']).write_text(json.dumps({\n"
        "    'argv': sys.argv[1:], 'cwd': os.getcwd(), 'cache': os.environ['XDG_CACHE_HOME'],\n"
        "    'job': os.environ.get('SLURM_JOB_ID')}))\n"
        f"sys.exit({exit_code})\n"
    )
    fake.chmod(0o755)
    env.update({"SNAKEMAKE_BIN": str(fake), "PHENORADAR_TEST_CAPTURE": str(capture)})
    if mode == "batch_cpu":
        del env["SLURM_MEM_PER_NODE"]
        env["SLURM_MEM_PER_CPU"] = "4096"
    arguments = ["--configfile", "config/data with spaces.yaml", "config/pilot.yaml",
                 "--config", "mem_gb=99",
                 "--", "prepare", "proteins"]
    cwd = checkout
    if mode == "direct":
        env = {k: v for k, v in env.items() if not k.startswith("SLURM_")}
        script = checkout / "run_pipeline.sh"
        # Ordinary execution must find the repository even from another directory.
        cwd = tmp_path
        arguments = ["--cores", "3", "--resources", "mem_gb=7", "disk_mb=12000",
                     "--software-deployment-method", "conda", *arguments]
    result = subprocess.run([str(script), *arguments], cwd=cwd,
                            env=env, capture_output=True, text=True)
    assert result.returncode == exit_code, result.stderr
    observed = json.loads(capture.read_text())
    argv = observed["argv"]
    assert argv[argv.index("--configfile") + 1:argv.index("--configfile") + 3] == [
        "config/data with spaces.yaml", "config/pilot.yaml"]
    assert argv[-3:] == ["--", "prepare", "proteins"]
    assert argv[argv.index("--executor") + 1] == "local"
    assert argv[argv.index("--software-deployment-method") + 1] == "conda"
    assert argv[argv.index("--config") + 1] == "mem_gb=99"
    assert argv[argv.index("--cores") + 1] == ("3" if mode == "direct" else "2")
    assert ("mem_mb=7000" if mode == "direct" else "mem_mb=4589") in argv
    if mode == "direct":
        assert "disk_mb=12000" in argv
        assert "mem_gb=7" not in argv
    else:
        assert "4.589 GB (4 GB reserved for overhead)" in result.stdout
    assert not any(arg.startswith("odb_slots=") for arg in argv)
    assert "--printshellcmds" in argv
    assert "--rerun-incomplete" in argv
    assert observed["cwd"] == str(checkout)
    assert observed["cache"] == str(checkout / ".cache")
    assert observed["job"] == (None if mode == "direct" else "123")


@pytest.mark.parametrize("changes,message", [
    ({"SLURM_JOB_NUM_NODES": "2"}, "Use one node and one task"),
    ({"SLURM_MEM_PER_NODE": "0"}, "--mem=0 is not supported"),
    ({"SLURM_MEM_PER_NODE": "2048"}, "Allocate more than 4 GB"),
])
def test_invalid_allocations_fail_before_work(batch_workspace, changes, message):
    checkout, script, env = batch_workspace
    env.update(changes)
    env["SNAKEMAKE_BIN"] = sys.executable
    result = subprocess.run(["bash", str(script)], cwd=checkout, env=env,
                            capture_output=True, text=True)
    assert result.returncode == 2
    assert message in result.stderr


@pytest.mark.parametrize("memory", ["0", "-1", "3.5", "128GB"])
def test_invalid_direct_memory_fails_before_work(batch_workspace, memory):
    checkout, _, env = batch_workspace
    env = {k: v for k, v in env.items() if not k.startswith("SLURM_")}
    env["SNAKEMAKE_BIN"] = sys.executable
    result = subprocess.run([str(checkout / "run_pipeline.sh"), "--resources", f"mem_gb={memory}"],
                            cwd=checkout, env=env, capture_output=True, text=True)
    assert result.returncode == 2
    assert "mem_gb must be a positive integer in GB" in result.stderr


@pytest.mark.parametrize("mode", ["batch", "direct", "direct_equals"])
def test_launcher_runs_dag_locally_without_submitting_jobs(batch_workspace, tmp_path, mode):
    checkout, script, env = batch_workspace
    snakemake = os.environ.get("SNAKEMAKE_BIN") or shutil.which("snakemake")
    if not snakemake:
        pytest.skip("Snakemake is not available")
    env["SNAKEMAKE_BIN"] = str(Path(snakemake).resolve())
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    marker = tmp_path / "unexpected-submission"
    for command in ("sbatch", "srun"):
        executable = fake_bin / command
        executable.write_text(f"#!{sys.executable}\nfrom pathlib import Path\n"
                              f"Path({str(marker)!r}).touch()\nraise SystemExit(99)\n")
        executable.chmod(0o755)
    env["PATH"] = str(fake_bin) + os.pathsep + env.get("PATH", "")
    (checkout / "workflow/Snakefile").write_text('''rule all:
    input: "a.txt", "b.txt"

rule task:
    output: "{sample}.txt"
    threads: 2
    resources: mem_mb=3000
    shell: "echo {threads} > {output}"
''')
    if mode == "batch":
        # Allocation settings must win over accidentally supplied larger budgets.
        arguments = ["--cores", "999", "--resources", "mem_gb=999"]
        cwd = checkout
        expected_cores, expected_memory = "2", "4589"
    else:
        env = {k: v for k, v in env.items() if not k.startswith("SLURM_")}
        script = checkout / "run_pipeline.sh"
        arguments = ["--cores", "1"]
        arguments += (["--resources=mem_gb=3"] if mode == "direct_equals"
                      else ["--resources", "mem_gb=3"])
        cwd = tmp_path
        expected_cores, expected_memory = "1", "3000"
    result = subprocess.run([str(script), *arguments], cwd=cwd, env=env,
                            capture_output=True, text=True, timeout=90)
    assert result.returncode == 0, result.stdout + result.stderr
    assert not marker.exists()
    assert (checkout / "a.txt").read_text().strip() == expected_cores
    assert (checkout / "b.txt").read_text().strip() == expected_cores
    assert f"Provided cores: {expected_cores}" in result.stdout + result.stderr
    assert f"mem_mb={expected_memory}" in result.stdout + result.stderr
    assert "odb_slots" not in result.stdout + result.stderr
