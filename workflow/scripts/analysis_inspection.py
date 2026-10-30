"""Read-only pending-work inspection through the pinned Snakemake DAG API."""
from collections import Counter
from dataclasses import replace
import contextlib
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
from unittest.mock import patch


def inspect_workflow(manifest, slurm, target):
    """Run the adapter in a child so Snakemake logging/cwd never affect the CLI."""
    root = Path(manifest['root'])
    environment = {**os.environ, 'PHENORADAR_READ_ONLY': '1', 'PYTHONDONTWRITEBYTECODE': '1'}
    environment.pop('PHENORADAR_FROZEN_IMAGE', None)
    if manifest.get('compatible_container_image'):
        environment['PHENORADAR_FROZEN_IMAGE'] = manifest['compatible_container_image']
    result = subprocess.run([sys.executable, str(Path(__file__).resolve()), '--inspect'],
        input=json.dumps({'root': str(root), 'analysis': manifest['pipeline']['output_root'],
                          'slurm': slurm, 'target': target,
                          'compatible_launcher': bool(manifest.get('compatible_launcher'))}),
        capture_output=True, text=True, timeout=300,
        env=environment)
    if result.returncode:
        raise ValueError('Cannot inspect pending work: ' + (result.stderr or result.stdout)[-3000:])
    return json.loads(result.stdout)


def plan(payload):
    # This is a deliberately narrow adapter to the host version pinned in
    # environment.yaml. It never instantiates an executor or runs a rule.
    import snakemake
    if snakemake.__version__ != '9.8.0':
        raise ValueError('Pending-work inspection requires Snakemake 9.8.0')
    from snakemake.api import SnakemakeApi
    from snakemake.common import async_run
    import snakemake.workflow as workflow_module
    from snakemake.persistence import Persistence
    from snakemake.settings.enums import Quietness, RerunTrigger
    from snakemake_interface_executor_plugins.settings import DeploymentMethod, SharedFSUsage
    from snakemake.settings.types import (ConfigSettings, DAGSettings, DeploymentSettings, ExecutionSettings,
        GroupSettings, OutputSettings, RemoteExecutionSettings, ResourceSettings, SchedulingSettings, StorageSettings)
    from snakemake.resources import DefaultResources
    from snakemake_interface_executor_plugins.registry import ExecutorPluginRegistry
    from alignment_batches import inspect_batches
    from phase_config import write_profile
    import yaml

    root = Path(payload['root']).resolve()
    analysis = root / payload['analysis']
    with tempfile.TemporaryDirectory(prefix='phenoradar-inspect-') as temporary:
        temporary = Path(temporary)
        profile = yaml.safe_load((write_profile(temporary / 'profile', payload['slurm']) / 'config.yaml').read_text())

        class InspectionPersistence(Persistence):
            def __init__(self, *args, **kwargs):
                # Writable caches/locks stay in temporary storage. Metadata and
                # incomplete records are read at their original paths/identities.
                super().__init__(*args, **dict(kwargs, path=temporary / 'persistence'))
                if (root / '.snakemake/metadata').is_dir():
                    self._metadata_path = str(root / '.snakemake/metadata')
                if (root / '.snakemake/incomplete').is_dir():
                    self._incomplete_path = str(root / '.snakemake/incomplete')
                # Environment identity checks must read the same deployment
                # paths as execution. A temporary prefix invents changed Conda
                # hashes and sends container queries to a nonexistent SIF.
                self.conda_env_path = str(root / '.snakemake/conda')
                self.conda_env_archive_path = str(root / '.snakemake/conda-archive')
                self.container_img_path = str(root / '.snakemake/singularity')

            def _record(self, *args, **kwargs):
                raise ValueError('Read-only inspection attempted to write workflow metadata')

            def _delete_record(self, *args, **kwargs):
                raise ValueError('Read-only inspection attempted to delete workflow metadata')

        resources = ResourceSettings(nodes=profile['jobs'],
            resources=profile.get('resources', {}),
            overwrite_threads=profile.get('set-threads', {}),
            overwrite_resources=profile.get('set-resources', {}),
            default_resources=DefaultResources(args=[f'{k}={v}' for k, v in profile['default-resources'].items()]))
        triggers = RerunTrigger.all()
        if payload.get('compatible_launcher'):
            triggers = triggers - {RerunTrigger.SOFTWARE_ENV}
        settings = DAGSettings(targets={payload['target']}, force_incomplete=True, rerun_triggers=triggers)
        with patch.object(workflow_module, 'Persistence', InspectionPersistence), \
                contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()), \
                SnakemakeApi(OutputSettings(quiet={Quietness.ALL})) as api, inspect_batches() as invalid:
            workflow_api = api.workflow(resource_settings=resources,
                config_settings=ConfigSettings(configfiles=[analysis / 'pipeline.yaml']),
                storage_settings=StorageSettings(shared_fs_usage=SharedFSUsage.all() - {SharedFSUsage.SOURCE_CACHE},
                                                 local_storage_prefix=temporary / 'storage'),
                deployment_settings=DeploymentSettings(
                    deployment_method={DeploymentMethod.CONDA, DeploymentMethod.APPTAINER}
                    if (root / 'VERSION').exists() else frozenset(), conda_base_path=Path('/opt/conda')),
                snakefile=root / 'workflow/Snakefile', workdir=root)
            workflow_api.dag(settings)
            workflow = workflow_api._workflow
            workflow.execution_settings = ExecutionSettings(lock=False)
            workflow.remote_execution_settings = RemoteExecutionSettings()
            workflow.scheduling_settings = SchedulingSettings(scheduler='greedy')
            workflow.group_settings = GroupSettings()
            plugin = ExecutorPluginRegistry().get_plugin('slurm')
            workflow._executor_plugin = replace(plugin, common_settings=replace(plugin.common_settings, dryrun_exec=True))
            for attempt in range(2):
                workflow._prepare_dag(forceall=False, ignore_incomplete=False, lock_warn_only=True, nolock=True)
                workflow._build_dag()
                async_run(workflow.dag.postprocess(update_needrun=False, check_initial=True))
                if not invalid or attempt:
                    break
                settings.forcerun = set(invalid)
            dag = workflow.dag
            pending, reused, reasons, examples = Counter(), Counter(), Counter(), []
            pending_jobs = list(dag.needrun_jobs())
            groups = set(dag._group.values())
            grouped = {job for group in groups for job in group.jobs}
            for job in pending_jobs:
                if not job.rule.norun:
                    dict(job.resources.items())  # Evaluate effective requests and budget guards.
            for group in groups:
                group.finalize()
                dict(group.resources.items())
            for job in dag.jobs:
                if job.rule.norun:
                    continue
                if not dag.needrun(job):
                    reused[job.rule.name] += 1
                    continue
                pending[job.rule.name] += 1
                reason = dag.reason(job)
                kinds = []
                if reason.missing_output:
                    kinds.append('missing_output')
                if reason.incomplete_output or any(async_run(workflow.persistence.incomplete(job))):
                    kinds.append('incomplete_output')
                if reason.params_changed or reason.code_changed or reason.software_stack_changed or reason.input_changed:
                    kinds.append('conditions_or_implementation_changed')
                if reason.updated_input:
                    kinds.append('updated_input')
                if any(str(p) in invalid for p in job.output):
                    kinds.append('invalid_batch_receipt')
                if not kinds:
                    kinds.append('dependency_or_forced')
                reasons.update(kinds)
                if len(examples) < 12:
                    examples.append({'rule': job.rule.name, 'reasons': kinds,
                                     'outputs': list(map(str, job.output))[:3]})
            unknown_checkpoints = [job.rule.name for job in pending_jobs if job.is_checkpoint]
            return {'pending': dict(pending), 'reused': dict(reused), 'reasons': dict(reasons),
                    'workflow_steps': sum(pending.values()),
                    'worker_jobs': len(groups) + sum(not j.rule.norun and not j.is_local and j not in grouped for j in pending_jobs),
                    'local_steps': sum(not j.rule.norun and j.is_local for j in pending_jobs),
                    'partial': bool(unknown_checkpoints), 'pending_checkpoints': unknown_checkpoints,
                    'examples': examples, 'invalid_batch_receipts': sorted(invalid)}


if __name__ == '__main__':
    try:
        print(json.dumps(plan(json.load(sys.stdin))))
    except Exception as error:
        print(str(error), file=sys.stderr)
        raise SystemExit(1)
