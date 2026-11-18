"""Charge effective per-job requests against the controller's remaining budget."""
import copy
from inspect import signature


BUDGETS = {'workflow_cpus': ('_cores', 'cpus', 1),
           'workflow_mem_mb': ('mem_mb', 'mem_gb', 1000)}


def apply_workflow_limits(workflow):
    """Attach global accounting resources after all rule overrides are applied.

    Separate accounting resources avoid Snakemake's automatic min(request, limit)
    adjustment of real CPU/memory requests. Validate before it can clip the charge.
    Resolve memory through a copy of the rule so dynamic resources keep Snakemake's
    normal input/wildcard/attempt semantics without recursive budget evaluation.
    """
    # Snakemake 9.8 accepts --local-cores without using it to schedule local
    # rules. An independent accounting resource enforces that limit without
    # capping the CPU requests or concurrency of remote workers.
    local_limit = workflow.resource_settings.local_cores
    if local_limit is not None:
        workflow.global_resources['controller_cpus'] = local_limit
        workflow.resource_scopes.update({'controller_cpus': 'global'})
        for rule in workflow.rules:
            if not rule.norun:
                rule.resources['controller_cpus'] = _local_charge(workflow, rule, local_limit)
    limits = {name: workflow.global_resources.get(name) for name in BUDGETS}
    limits = {name: value for name, value in limits.items() if value is not None}
    if not limits: return
    workflow.resource_scopes.update({name: 'global' for name in limits})
    for rule in workflow.rules:
        # Dependency-only targets allocate no worker resources.
        if rule.norun: continue
        probe = copy.copy(rule)
        probe.resources = {key: value for key, value in rule.resources.items()
                           if key in {'_cores', 'mem_mb', 'mem'}}
        for name, limit in limits.items():
            rule.resources[name] = _charge(probe, name, limit)


def _local_charge(workflow, rule, limit):
    def charge(wildcards, threads):
        # The executor is selected after the Snakefile is parsed.
        if workflow.local_exec or not workflow.is_local(rule):
            return 0
        if threads > limit:
            raise ValueError(f'{rule.name} requires {threads} CPUs, exceeding the controller CPU limit of {limit}')
        return threads
    return charge


def _charge(probe, name, limit):
    resource, setting, scale = BUDGETS[name]

    def checked(value):
        if type(value) is not int or value < 1:
            raise ValueError(f'{probe.name} must declare a positive {setting} request to use total_limits.{setting}')
        if value > limit:
            raise ValueError(f'{probe.name} requires {value / scale:g} {setting}, exceeding the remaining '
                             f'total_limits.{setting} budget of {limit / scale:g} after the controller allocation')
        return value

    if resource == '_cores':
        return lambda wildcards, threads: checked(threads)
    value = probe.resources.get(resource)
    if type(value) is int:
        return lambda wildcards: checked(value)
    if callable(value) and tuple(signature(value).parameters) == ('wildcards',):
        # Wildcard-only requests can be charged when a group is first planned,
        # even if its intermediate input files have not been produced yet.
        return lambda wildcards, attempt: checked(probe.expand_resources(wildcards, None, attempt).get(resource))

    def charge(wildcards, input, threads, attempt):
        return checked(probe.expand_resources(wildcards, input, attempt).get(resource))

    return charge
