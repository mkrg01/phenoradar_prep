"""Use the installed Snakemake methods that failed before BUSCO execution."""
import asyncio
from types import SimpleNamespace

import pytest

from slurm_submit import configure_group_wait_files


@pytest.fixture
def grouped_inputs(tmp_path, monkeypatch):
    pytest.importorskip('snakemake')
    import snakemake.api
    from snakemake.jobs import GroupJob, Job
    from snakemake.io import IOFile
    from snakemake_interface_executor_plugins.settings import DeploymentMethod, ExecMode
    # Restore the actual class method after each test, including installation.
    monkeypatch.setattr(GroupJob, 'get_wait_for_files', GroupJob.get_wait_for_files)
    monkeypatch.setattr(snakemake, '__version__', '9.8.0')
    workflow = SimpleNamespace(
        exec_mode=ExecMode.DEFAULT, executor_plugin=SimpleNamespace(name='slurm'),
        deployment_settings=SimpleNamespace(deployment_method={DeploymentMethod.CONDA,
                                                               DeploymentMethod.APPTAINER}))
    dag = SimpleNamespace(workflow=workflow)
    source = tmp_path / 'marker.faa'
    source.write_text('>sample\nMPEPTIDE\n')
    intermediate = IOFile(str(tmp_path / 'alignment.faa'))
    env = SimpleNamespace(address=str(tmp_path / 'container-only/conda-env'),
                          is_containerized=True, is_externally_managed=False)
    def job(inputs):
        return SimpleNamespace(local_input=inputs, storage_input=[], shadow_dir=None,
                               conda_env=env, dag=dag)
    group = SimpleNamespace(jobs=[job([IOFile(str(source))]), job([intermediate])],
                            all_products={intermediate}, dag=dag)
    return workflow, group, env, GroupJob, Job


def test_containerized_group_waits_succeed_and_internal_products_stay_excluded(grouped_inputs):
    from snakemake.io import wait_for_files
    workflow, group, env, GroupJob, Job = grouped_inputs
    upstream = GroupJob.get_wait_for_files(group)
    assert env.address in upstream
    assert all(env.address not in Job.get_wait_for_files(job) for job in group.jobs)
    assert configure_group_wait_files(workflow)
    wait = GroupJob.get_wait_for_files(group)
    assert env.address not in wait and not set(wait) & group.all_products
    asyncio.run(wait_for_files(wait, latency_wait=0))


@pytest.mark.parametrize('containerized,managed,expected', [
    (False, False, True), (False, True, False), (True, True, False)])
def test_host_and_externally_managed_environments_keep_single_job_semantics(
        grouped_inputs, containerized, managed, expected):
    workflow, group, env, GroupJob, _ = grouped_inputs
    env.is_containerized, env.is_externally_managed = containerized, managed
    configure_group_wait_files(workflow)
    assert (env.address in GroupJob.get_wait_for_files(group)) == expected


def test_missing_inputs_and_shadow_directories_are_still_required(grouped_inputs):
    workflow, group, env, GroupJob, _ = grouped_inputs
    # Even a real input with the same spelling as a container path must remain.
    group.jobs[0].local_input.append(env.address)
    group.jobs[0].shadow_dir = env.address + '-shadow'
    configure_group_wait_files(workflow)
    wait = GroupJob.get_wait_for_files(group)
    assert env.address in wait and group.jobs[0].shadow_dir in wait


def test_other_workflows_are_unchanged_and_installation_is_idempotent(grouped_inputs):
    workflow, group, env, GroupJob, _ = grouped_inputs
    configure_group_wait_files(workflow)
    installed = GroupJob.get_wait_for_files
    configure_group_wait_files(workflow)
    assert GroupJob.get_wait_for_files is installed
    group.dag = SimpleNamespace(workflow=SimpleNamespace(deployment_settings=workflow.deployment_settings))
    assert env.address in GroupJob.get_wait_for_files(group)


@pytest.mark.parametrize('change', ['version', 'remote', 'local', 'no_container', 'no_conda'])
def test_compatibility_is_limited_to_pinned_containerized_slurm_controller(grouped_inputs, monkeypatch, change):
    import snakemake
    from snakemake_interface_executor_plugins.settings import DeploymentMethod, ExecMode
    workflow, _, _, GroupJob, _ = grouped_inputs
    before = GroupJob.get_wait_for_files
    if change == 'version': monkeypatch.setattr(snakemake, '__version__', '9.9.0')
    elif change == 'remote': workflow.exec_mode = ExecMode.REMOTE
    elif change == 'local': workflow.executor_plugin.name = 'local'
    elif change == 'no_container': workflow.deployment_settings.deployment_method.remove(DeploymentMethod.APPTAINER)
    elif change == 'no_conda': workflow.deployment_settings.deployment_method.remove(DeploymentMethod.CONDA)
    assert not configure_group_wait_files(workflow)
    assert GroupJob.get_wait_for_files is before
