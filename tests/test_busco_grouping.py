"""Exercise BUSCO grouping and resource accounting with Snakemake 9.8."""
from pathlib import Path
from unittest.mock import patch

import pytest
import yaml

from phase_config import merge_slurm, write_profile

ROOT = Path(__file__).resolve().parents[1]
BRANCHES = ("phylogeny/all", "phylogeny/representatives")
MARKERS = ("1at1", "2at1")
STAGES = ("align_busco_marker", "trim_busco_marker", "infer_busco_gene_tree")


@pytest.fixture
def marker_workflow(tmp_path):
    # Use the production rules, with already-collected tiny marker inputs.
    # No external alignment tools, containers, or Slurm commands are executed.
    rules = (ROOT / "workflow/rules/phylogeny.smk").read_text()
    rules = "rule align_busco_marker:" + rules.split("rule align_busco_marker:", 1)[1]
    rules = rules.split("rule merge_busco_gene_trees:", 1)[0]
    prefix = f"""
import json
import sys
sys.path.insert(0, {str(ROOT / "workflow/scripts")!r})
OUT = "results/test"
LOG = "logs/test"
SCRIPTS = {str(ROOT / "workflow/scripts")!r}
PHYLO_RUN = OUT + "/{{phylo_branch}}"
PHY = {{"min_taxa": 4, "min_protein_length": 100, "trimal_mode": "gappyout"}}
PYTHON = "python"
RANDOM_SEED = 12345

rule all:
    input:
        expand(OUT + "/{{phylo_branch}}/gene_trees/{{marker}}.nwk",
               phylo_branch={BRANCHES!r}, marker={MARKERS!r})

rule collect_busco_markers:
    output: fasta=directory(PHYLO_RUN + "/markers")
"""
    suffix = """
from resource_limits import apply_workflow_limits
apply_workflow_limits(workflow)
"""
    (tmp_path / "Snakefile").write_text(prefix + rules + suffix)
    for branch in BRANCHES:
        out = tmp_path / "results/test" / branch
        (out / "plan").mkdir(parents=True)
        (out / "plan/species.tsv").write_text("species\nPlant_A\n")
        (out / "markers").mkdir()
        for marker in MARKERS:
            (out / "markers" / (marker + ".faa")).write_text(">Plant_A_g1\nMPEPTIDE\n")
        (out / "markers/.snakemake_timestamp").touch()
    return tmp_path


def inspect_groups(project, overrides=None):
    pytest.importorskip("snakemake")
    from snakemake.api import SnakemakeApi
    from snakemake.executors import dryrun
    from snakemake.resources import DefaultResources
    from snakemake.settings.types import DAGSettings, ResourceSettings

    public = yaml.safe_load((ROOT / "config/analysis.yaml").read_text())["slurm"]
    slurm = merge_slurm(public, overrides or {}, build=False)
    profile_dir = write_profile(project / "profile", slurm)
    profile = yaml.safe_load((profile_dir / "config.yaml").read_text())
    resources = ResourceSettings(
        cores=64, nodes=profile["jobs"],
        resources=profile["resources"],
        overwrite_threads=profile["set-threads"],
        overwrite_resources=profile["set-resources"],
        default_resources=DefaultResources(
            args=[f"{key}={value}" for key, value in profile["default-resources"].items()]),
    )
    # Use the dry-run executor with cluster resource aggregation. The default
    # local dry-run mode sums resources across sequential group blocks.
    with patch.object(dryrun.common_settings, "non_local_exec", True), SnakemakeApi() as api:
        workflow = api.workflow(resource_settings=resources, snakefile=project / "Snakefile",
                                workdir=project)
        targets = {f"results/test/{branch}/gene_trees/{marker}.nwk"
                   for branch in BRANCHES for marker in MARKERS}
        workflow.dag(DAGSettings(targets=targets)).execute_workflow(executor="dryrun")
        groups = set(workflow._workflow.dag._group.values())
        return [{
            "members": {(job.wildcards.phylo_branch, job.wildcards.marker) for job in group.jobs},
            "rules": {job.rule.name for job in group.jobs},
            "resources": dict(group.resources.items()),
        } for group in groups]


def test_one_group_per_marker_and_tree_set_with_serial_resource_limits(marker_workflow):
    # Exactly one default group fits: series must take max(CPU/memory), not sum.
    groups = inspect_groups(marker_workflow, {"total_limits": {"cpus": 5, "mem_gb": 16}})
    assert len(groups) == len(BRANCHES) * len(MARKERS)
    assert {next(iter(group["members"])) for group in groups} == {
        (branch, marker) for branch in BRANCHES for marker in MARKERS}
    for group in groups:
        assert len(group["members"]) == 1
        assert group["rules"] == set(STAGES)
        resources = group["resources"]
        assert resources["_cores"] == resources["workflow_cpus"] == 4
        assert resources["mem_mb"] == resources["workflow_mem_mb"] == 8000
        assert resources["runtime"] == 3 * 24 * 60


def test_group_honors_per_stage_resource_overrides(marker_workflow):
    groups = inspect_groups(marker_workflow, {"per_job_resources": {
        STAGES[0]: {"cpus": 6, "mem_gb": 12, "time": "01:30:00"},
        STAGES[1]: {"time": "00:30:00"},
        STAGES[2]: {"cpus": 8, "mem_gb": 16, "time": "02:00:00"},
    }})
    for group in groups:
        resources = group["resources"]
        assert resources["_cores"] == resources["workflow_cpus"] == 8
        assert resources["mem_mb"] == resources["workflow_mem_mb"] == 16000
        assert resources["runtime"] == 240


@pytest.mark.parametrize("completed", [1, 2, 3])
def test_resume_only_groups_unfinished_stages(marker_workflow, completed):
    out = marker_workflow / "results/test" / BRANCHES[0]
    outputs = [
        [out / "alignments/raw/1at1.faa", out / "alignments/raw/1at1.json"],
        [out / "alignments/1at1.faa", out / "alignments/1at1.json",
         out / "alignments/1at1.columns.tsv"],
        [out / "gene_trees/1at1.nwk", out / "gene_trees/1at1.json"],
    ]
    saved = {}
    for stage in outputs[:completed]:
        for path in stage:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("completed stage\n")
            saved[path] = (path.read_bytes(), path.stat().st_mtime_ns)
    groups = inspect_groups(marker_workflow)
    selected = [group for group in groups if (BRANCHES[0], MARKERS[0]) in group["members"]]
    if completed == 3:
        assert not selected
    else:
        assert len(selected) == 1
        assert selected[0]["rules"] == set(STAGES[completed:])
        assert selected[0]["resources"]["runtime"] == (3 - completed) * 24 * 60
    assert all((path.read_bytes(), path.stat().st_mtime_ns) == state for path, state in saved.items())
