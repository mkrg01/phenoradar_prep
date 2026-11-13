"""Public build/analysis settings and per-submission Slurm resources."""
import copy
import os
import re
import subprocess
import shlex
import sys
from pathlib import Path
import yaml


def read_yaml(path):
    value = yaml.safe_load(Path(path).read_text())
    if not isinstance(value, dict):
        raise ValueError(f"configuration must be a mapping: {path}")
    return value


DEFAULT_REGISTERED_SUPPORT = 0.10
DEFAULT_TAXONOMY_BATCH_SIZE = 24
TAXONOMY_SCREEN_METHOD = 'angiosperms353'


def taxonomy_screen_settings(value):
    defaults = {'enabled': False, 'method': TAXONOMY_SCREEN_METHOD,
                'max_registered_support': DEFAULT_REGISTERED_SUPPORT,
                'batch_size': DEFAULT_TAXONOMY_BATCH_SIZE}
    if not isinstance(value, dict) or set(value) - defaults.keys():
        raise ValueError('taxonomy_screen accepts enabled, method, max_registered_support and batch_size')
    result = {**defaults, **value}
    if type(result['enabled']) is not bool:
        raise ValueError('taxonomy_screen.enabled must be true or false')
    if result['method'] != TAXONOMY_SCREEN_METHOD:
        raise ValueError(f'taxonomy_screen.method must be {TAXONOMY_SCREEN_METHOD}')
    cutoff = result['max_registered_support']
    if type(cutoff) not in (int, float) or not 0 <= cutoff <= 1:
        raise ValueError('taxonomy_screen.max_registered_support must be between zero and one')
    if type(result['batch_size']) is not int or result['batch_size'] < 1:
        raise ValueError('taxonomy_screen.batch_size must be a positive integer')
    return result


def deep_merge(base, override):
    result = copy.deepcopy(base)
    for key, value in override.items():
        result[key] = deep_merge(result[key], value) if isinstance(value, dict) and isinstance(result.get(key), dict) else copy.deepcopy(value)
    return result


def time_minutes(value):
    """Convert a quoted Slurm duration to minutes, rounding seconds up."""
    error = f'invalid Slurm time: {value!r}; use a positive quoted duration such as "1-00:00:00"'
    if not isinstance(value, str) or not re.fullmatch(r'(?:[0-9]+-)?[0-9]+(?::[0-9]+){0,2}', value):
        raise ValueError(error)
    day_text, separator, clock = value.partition('-')
    days = int(day_text) if separator else 0
    parts = [int(part) for part in (clock if separator else day_text).split(':')]
    if separator or len(parts) == 3:
        hours, minutes, seconds = parts + [0] * (3 - len(parts))
    else:
        hours = 0
        minutes = parts[0]
        seconds = parts[1] if len(parts) == 2 else 0
    total_seconds = days * 86400 + hours * 3600 + minutes * 60 + seconds
    if total_seconds <= 0:
        raise ValueError(error)
    return (total_seconds + 59) // 60


# Normalize the public layout to the historical execution schema. Saved manifests
# and execution receipts retain their original shape and checksums on disk.
DIRECT_JOBS = {'sample', 'assembly', 'busco', 'quant', 'controller', 'taxonomy_screen'}
WORKER_DEFAULTS = {'mem_gb': 8, 'time': '1-00:00:00'}
BUDGET_RESOURCES = {'cpus': 'workflow_cpus', 'mem_gb': 'workflow_mem_mb'}


def execution_config(config):
    """The public execution section replaces slurm; old saved settings still load."""
    if 'execution' in config and 'slurm' in config:
        raise ValueError('use execution or legacy slurm settings, not both')
    return config.get('execution', config.get('slurm', {}))


def normalize_slurm(slurm, build=True):
    if not isinstance(slurm, dict): raise ValueError('slurm must be a mapping')
    result = copy.deepcopy(slurm)
    # Old build snapshots/retry files may contain this retired setting. Array
    # limits now always come from Slurm; leave the saved inputs unchanged.
    if build: result.pop('array_size', None)
    if 'total_limits' in result:
        limits = result['total_limits']
        if not isinstance(limits, dict) or set(limits) - {'jobs', 'cpus', 'mem_gb'}:
            raise ValueError('invalid slurm.total_limits')
        if 'jobs' in limits:
            if 'jobs' in result or 'concurrency' in result:
                raise ValueError('do not mix total_limits.jobs with legacy jobs/concurrency')
            result['jobs'] = limits.pop('jobs')
            if build: result['concurrency'] = result['jobs']
    if 'per_job_resources' in result:
        resources = result.pop('per_job_resources')
        if not isinstance(resources, dict): raise ValueError('slurm.per_job_resources must be a mapping')
        for name, values in resources.items():
            section = 'stages' if name in DIRECT_JOBS else 'rules'
            target = result.setdefault(section, {})
            if not isinstance(target, dict): raise ValueError(f'slurm.{section} must be a mapping')
            if name in target: raise ValueError(f'duplicate resources for {name}')
            target[name] = values
    if build and isinstance(result.get('stages'), dict):
        stages = result['stages']
        legacy = set(stages) & {'assembly', 'busco', 'quant'}
        if legacy:
            if 'sample' in stages:
                raise ValueError('do not mix sample and legacy assembly/busco/quant resources')
            if legacy == {'assembly', 'busco', 'quant'}:
                for name in legacy:
                    job = stages[name]
                    if not isinstance(job, dict) or set(job) != {'cpus', 'mem_gb', 'time'}:
                        raise ValueError(f'invalid resources for {name}')
                    for key in ('cpus', 'mem_gb'):
                        if type(job[key]) is not int or job[key] < 1:
                            raise ValueError(f'invalid {name}.{key}')
                minutes = sum(time_minutes(stages[name]['time']) for name in legacy)
                days, remainder = divmod(minutes, 1440)
                hours, minute = divmod(remainder, 60)
                stages['sample'] = {'cpus': max(stages[name]['cpus'] for name in legacy),
                                    'mem_gb': max(stages[name]['mem_gb'] for name in legacy),
                                    'time': f'{days}-{hours:02d}:{minute:02d}:00'}
                for name in legacy: del stages[name]
            elif legacy == {'assembly'}:
                # Existing assembly-only retry overrides now size the whole sample job.
                stages['sample'] = stages.pop('assembly')
            else:
                raise ValueError('use per_job_resources.sample for full-sample retry resources')
    return result


def merge_slurm(base, override, build=True):
    # Normalize each side before merging so old retry files can override new configs
    # and new retry files can override old snapshots (including unequal old caps).
    base = normalize_slurm(base, build)
    base.setdefault('default_resources', copy.deepcopy(WORKER_DEFAULTS))
    return validate_slurm(deep_merge(base, normalize_slurm(override, build)), build)


def validate_slurm(slurm, build=True):
    slurm = normalize_slurm(slurm, build)
    slurm.setdefault('default_resources', copy.deepcopy(WORKER_DEFAULTS))
    slurm.setdefault('rules', {})
    if slurm.get('partition') is not None and (not isinstance(slurm['partition'], str) or not slurm['partition'].strip()):
        raise ValueError('slurm.partition must be null or a nonempty string')
    allowed = {'partition','jobs','stages','default_resources','rules','total_limits','scheduler','age'}
    if build: allowed.add('concurrency')
    if set(slurm) - allowed: raise ValueError('unknown slurm setting')
    if slurm.get('scheduler', 'auto') not in {'auto', 'slurm', 'age'}:
        raise ValueError('execution.scheduler must be auto, slurm or age')
    age = slurm.get('age', {})
    if not isinstance(age, dict) or set(age) - {'parallel_environment', 'memory_resource', 'queue', 'extra_resources'}:
        raise ValueError('unknown execution.age settings')
    for key in ('parallel_environment', 'memory_resource', 'queue'):
        if age.get(key) is not None and (not isinstance(age[key], str) or not re.fullmatch(r'[A-Za-z0-9_.@*,-]+', age[key])):
            raise ValueError(f'invalid execution.age.{key}')
    extra = age.get('extra_resources', {})
    if not isinstance(extra, dict) or any(not re.fullmatch(r'[A-Za-z][A-Za-z0-9_]*', k) or
            not (v is True or isinstance(v, (int, float, str)) and re.fullmatch(r'[A-Za-z0-9_.+-]+', str(v))) for k, v in extra.items()):
        raise ValueError('invalid execution.age.extra_resources')
    for key in ('jobs', 'concurrency') if build else ('jobs',):
        if key not in slurm or (slurm[key] is not None and (type(slurm[key]) is not int or slurm[key] < 1)):
            raise ValueError(f'slurm.{key} must be a positive integer or null')
    for key, value in slurm.get('total_limits', {}).items():
        if value is not None and (type(value) is not int or value < 1):
            raise ValueError(f'slurm.total_limits.{key} must be a positive integer or null')
    expected = {'sample','controller'} if build else {'controller'}
    allowed_stages = expected | {'taxonomy_screen'} if build else expected
    if not isinstance(slurm.get('stages'), dict) or not expected <= slurm['stages'].keys() <= allowed_stages:
        raise ValueError(f'slurm.per_job_resources must contain {sorted(expected)}; allowed stages: {sorted(allowed_stages)}')
    for stage, job in slurm['stages'].items():
        if not isinstance(job, dict) or set(job) != {'cpus','mem_gb','time'}:
            raise ValueError(f'invalid resources for {stage}')
        for key in ('cpus','mem_gb'):
            if type(job[key]) is not int or job[key] < 1: raise ValueError(f'invalid {stage}.{key}')
        time_minutes(job['time'])
    if not isinstance(slurm.get('rules', {}), dict): raise ValueError('slurm.rules must be a mapping')
    if set(slurm.get('rules', {})) & DIRECT_JOBS: raise ValueError('direct jobs cannot be Snakemake rule overrides')
    for section, values in [('default_resources', slurm['default_resources']), *slurm.get('rules', {}).items()]:
        keys = {'mem_gb','time'} if section == 'default_resources' else {'cpus','mem_gb','time'}
        if not isinstance(values, dict) or set(values) - keys: raise ValueError(f'invalid rule resources: {section}')
        for key, value in values.items():
            if key == 'time': time_minutes(value)
            elif type(value) is not int or value < 1: raise ValueError(f'invalid {section}.{key}')
    return slurm


def execution_settings(frozen, resources=None, build=True):
    # Biological settings are never read from a resubmission resource override.
    override = read_yaml(resources) if resources else {}
    if override and not ({'slurm', 'execution'} & override.keys()):
        raise ValueError('resource override requires an execution or legacy slurm section')
    return merge_slurm(frozen, execution_config(override), build)


def resolve_array_size():
    """Read the current Slurm limits for one-based sample arrays."""
    try:
        output = subprocess.check_output(['scontrol', 'show', 'config'], text=True,
                                         stderr=subprocess.PIPE, timeout=10)
    except (OSError, subprocess.SubprocessError) as error:
        raise ValueError("Cannot query Slurm array limits with 'scontrol show config'; "
                         "check Slurm availability and retry") from error
    settings = {}
    for line in output.splitlines():
        name, separator, raw = line.partition('=')
        if separator:
            settings[name.strip()] = raw.strip()
    raw = settings.get('MaxArraySize', '')
    if not re.fullmatch(r'[0-9]+', raw) or 'SchedulerParameters' not in settings:
        raise ValueError('Cannot read Slurm MaxArraySize/SchedulerParameters from scontrol show config')
    # Slurm's maximum task ID is MaxArraySize - 1; our workers start at one.
    size = int(raw) - 1
    if size < 1:
        raise ValueError('Slurm MaxArraySize must be at least 2 for one-based sample arrays')
    for parameter in settings['SchedulerParameters'].split(','):
        key, separator, raw = parameter.strip().partition('=')
        if key != 'max_array_tasks':
            continue
        if not separator or not re.fullmatch(r'[0-9]+', raw) or int(raw) < 1:
            raise ValueError('Cannot read a positive Slurm max_array_tasks limit')
        size = min(size, int(raw))
    return size


def array_concurrency(slurm, stage):
    """Limit running array tasks by job count and the per-run request budgets."""
    result = slurm['concurrency']
    job = slurm['stages'][stage]
    for key, total in slurm.get('total_limits', {}).items():
        if total is not None:
            cap = total // job[key]
            result = cap if result is None else min(result, cap)
            if total < job[key]:
                raise ValueError(f'{stage} requires {job[key]} {key}, exceeding total_limits.{key}={total}')
    return result


def workflow_jobs(slurm):
    """Translate the public null limit to Snakemake's explicit unlimited value."""
    return 'unlimited' if slurm['jobs'] is None else slurm['jobs']


def worker_budgets(slurm):
    """Reserve the controller allocation before scheduling worker jobs."""
    result = {}
    for key, total in slurm.get('total_limits', {}).items():
        if total is None: continue
        remaining = total - slurm['stages']['controller'][key]
        if remaining < 1:
            raise ValueError(f'total_limits.{key}={total} leaves no worker budget after the controller allocation')
        result[BUDGET_RESOURCES[key]] = remaining * (1000 if key == 'mem_gb' else 1)
    return result


def _snakemake_resources(values):
    resources = {key:value for key,value in values.items() if key != 'cpus'}
    if 'mem_gb' in resources:
        # Preserve the existing Slurm requests in MB.
        resources['mem_mb'] = resources.pop('mem_gb') * 1000
    if 'time' in resources:
        resources['runtime'] = time_minutes(resources.pop('time'))
    return resources


def write_profile(destination, slurm, *, root=None, pipeline=None, usage=None, inspection=False):
    # Also used with partial legacy resource dictionaries by low-level callers.
    # Snakemake profiles only use worker limits, not sample-array concurrency.
    slurm = normalize_slurm(slurm, build=False)
    defaults = _snakemake_resources(slurm.get('default_resources', WORKER_DEFAULTS))
    if slurm.get('partition'): defaults['slurm_partition'] = slurm['partition']
    profile = {'executor':'slurm', 'jobs':workflow_jobs(slurm), 'latency-wait':90,
               'rerun-incomplete':True, 'printshellcmds':True, 'default-resources':defaults,
               'slurm-no-account':True,
               'slurm-status-command':'sacct', 'slurm-init-seconds-before-status-checks':10,
               'seconds-between-status-checks':10}
    # run_pipeline.sh installs the serialized sbatch gate. Its two-second limit
    # does not need Snakemake's accumulating scheduler reevaluation timers.
    profile['scheduler'] = 'greedy'
    profile['set-resources'] = {rule:_snakemake_resources(values) for rule,values in slurm.get('rules', {}).items()}
    profile['set-threads'] = {rule:values['cpus'] for rule,values in slurm.get('rules', {}).items() if 'cpus' in values}
    profile['set-resources'] = {rule: values for rule, values in profile['set-resources'].items() if values}
    for key in ('set-resources', 'set-threads'):
        if not profile[key]: del profile[key]
    budgets = worker_budgets(slurm)
    if budgets:
        profile['resources'] = budgets
        profile['set-resource-scopes'] = {key: 'global' for key in budgets}
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    if not inspection and slurm.get('scheduler') == 'age':
        from common import write_json
        root = Path(root or os.getcwd()).resolve()
        context = {'execution': slurm, 'root': str(root),
                   'pipeline': str(pipeline or root / 'workflow/pipeline_defaults.yaml'),
                   'usage': str(usage or root / '.snakemake/storage')}
        settings = destination / 'execution.json'
        write_json(settings, context)
        script = Path(__file__).with_name('scheduler.py')
        profile = {key: value for key, value in profile.items() if not key.startswith('slurm-')}
        profile.update(executor='cluster-generic',
                       **{'cluster-generic-submit-cmd': shlex.join([sys.executable, str(script), 'submit-job', '--settings', str(settings)]),
                          'cluster-generic-status-cmd': shlex.join([sys.executable, str(script), 'status-job', '--state', str(root / '.snakemake/cluster_submissions')]),
                          'cluster-generic-cancel-cmd': 'qdel'})
    (destination / 'config.yaml').write_text(yaml.safe_dump(profile, sort_keys=False))
    return destination
