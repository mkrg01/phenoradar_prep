"""Public build/analysis settings and per-submission Slurm resources."""
import copy
import re
from pathlib import Path
import yaml


def read_yaml(path):
    value = yaml.safe_load(Path(path).read_text())
    if not isinstance(value, dict):
        raise ValueError(f"configuration must be a mapping: {path}")
    return value


def deep_merge(base, override):
    result = copy.deepcopy(base)
    for key, value in override.items():
        result[key] = deep_merge(result[key], value) if isinstance(value, dict) and isinstance(result.get(key), dict) else copy.deepcopy(value)
    return result


def validate_slurm(slurm, build=True):
    if not isinstance(slurm, dict): raise ValueError('slurm must be a mapping')
    for key in ('partition', 'account'):
        if slurm.get(key) is not None and (not isinstance(slurm[key], str) or not slurm[key].strip()):
            raise ValueError(f'slurm.{key} must be null or a nonempty string')
    allowed = {'partition','account','jobs','stages','default_resources','rules'}
    if build: allowed.update({'concurrency','array_size'})
    if set(slurm) - allowed: raise ValueError('unknown slurm setting')
    for key in ('jobs', 'concurrency', 'array_size') if build else ('jobs',):
        if type(slurm.get(key)) is not int or slurm[key] < 1: raise ValueError(f'slurm.{key} must be positive')
    expected = {'assembly','busco','quant','controller'} if build else {'controller'}
    if set(slurm.get('stages', {})) != expected: raise ValueError(f'slurm.stages must contain {sorted(expected)}')
    for stage, job in slurm['stages'].items():
        if set(job) != {'cpus','mem_mb','time'}: raise ValueError(f'invalid resources for {stage}')
        for key in ('cpus','mem_mb'):
            if type(job[key]) is not int or job[key] < 1: raise ValueError(f'invalid {stage}.{key}')
        if not re.fullmatch(r'(?:[0-9]+-)?[0-9]+(?::[0-9]{2}){0,2}', str(job['time'])): raise ValueError(f'invalid Slurm time: {stage}')
    for section, values in [('default_resources', slurm.get('default_resources', {})), *slurm.get('rules', {}).items()]:
        keys = {'mem_mb','runtime'} if section == 'default_resources' else {'cpus','mem_mb','runtime'}
        if not isinstance(values, dict) or set(values) - keys: raise ValueError(f'invalid rule resources: {section}')
        for key, value in values.items():
            if type(value) is not int or value < 1: raise ValueError(f'invalid {section}.{key}')
    return slurm


def execution_settings(frozen, resources=None, build=True):
    # Biological settings are never read from a resubmission resource override.
    override = read_yaml(resources) if resources else {}
    if override and 'slurm' not in override: raise ValueError('resource override requires a slurm section')
    return validate_slurm(deep_merge(frozen, override.get('slurm', {})), build)


def write_profile(destination, slurm):
    defaults = dict(slurm.get('default_resources', {}))
    defaults.pop('cpus', None)
    for key in ('partition','account'):
        if slurm.get(key): defaults['slurm_' + key] = slurm[key]
    profile = {'executor':'slurm', 'jobs':slurm['jobs'], 'latency-wait':90,
               'rerun-incomplete':True, 'printshellcmds':True, 'default-resources':defaults}
    profile['set-resources'] = {rule:{k:v for k,v in values.items() if k != 'cpus'} for rule,values in slurm.get('rules', {}).items()}
    profile['set-threads'] = {rule:values['cpus'] for rule,values in slurm.get('rules', {}).items() if 'cpus' in values}
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    (destination / 'config.yaml').write_text(yaml.safe_dump(profile, sort_keys=False))
    return destination
