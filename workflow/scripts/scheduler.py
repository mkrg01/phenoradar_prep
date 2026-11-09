"""Scheduler detection and shared Slurm/AGE submission and accounting adapters."""
import argparse
import copy
import hashlib
from functools import lru_cache
import json
import math
import os
from pathlib import Path
import pwd
import re
import shlex
import shutil
import subprocess
import sys
import time
import uuid
import xml.etree.ElementTree as ET

from common import now, write_json


def detect(requested='auto', *, required=True):
    requested = os.environ.get('PHENORADAR_SCHEDULER') or requested
    if requested not in {'auto', 'slurm', 'age'}:
        raise ValueError('scheduler must be auto, slurm or age')
    if requested != 'auto':
        return requested
    allocated = []
    if os.environ.get('SLURM_JOB_ID'):
        allocated.append('slurm')
    if os.environ.get('JOB_ID') and (os.environ.get('SGE_ROOT') or os.environ.get('NSLOTS')):
        allocated.append('age')
    candidates = allocated or [name for name, commands in (
        ('slurm', ('sbatch', 'squeue', 'sacct')), ('age', ('qsub', 'qstat', 'qacct')))
        if all(shutil.which(command) for command in commands)]
    if len(candidates) == 1:
        return candidates[0]
    if required:
        reason = 'Both Slurm and AGE are available' if candidates else 'No scheduler client found'
        raise ValueError(reason + '; set execution.scheduler or PHENORADAR_SCHEDULER explicitly')
    return None


def current_id():
    return os.environ.get('SLURM_JOB_ID') or os.environ.get('JOB_ID')


@lru_cache(maxsize=1)
def age_capabilities():
    """Read the site's PE and memory complexes rather than assume a host name."""
    try:
        pes = subprocess.check_output(['qconf', '-spl'], text=True, timeout=10).split()
        complexes = subprocess.check_output(['qconf', '-sc'], text=True, timeout=10)
    except (OSError, subprocess.SubprocessError) as error:
        raise ValueError('Cannot inspect AGE resources with qconf; set execution.age parallel_environment and memory_resource') from error
    names = {line.split()[0] for line in complexes.splitlines() if line.strip() and not line.startswith('#')}
    pe = next((p for p in ('def_slot', 'smp', 'threads') if p in pes), pes[0] if len(pes) == 1 else None)
    memory = next((m for m in ('s_vmem', 'h_vmem', 'mem_req') if m in names), None)
    if not pe or not memory:
        raise ValueError('AGE resources are ambiguous; set execution.age parallel_environment and memory_resource')
    # Every tool in this workflow is a shared-memory process on one node.
    description = subprocess.check_output(['qconf', '-sp', pe], text=True, timeout=10)
    if not re.search(r'^allocation_rule\s+\$pe_slots\s*$', description, re.M):
        raise ValueError(f'AGE PE {pe} is not confined to one node; configure a shared-memory parallel_environment')
    return pe, memory


def allocation():
    """Read an AGE allocation for low-level run_pipeline.sh use inside qsub."""
    slots = int(os.environ.get('NSLOTS', '1'))
    if os.environ.get('PE_HOSTFILE'):
        hosts = {line.split()[0] for line in Path(os.environ['PE_HOSTFILE']).read_text().splitlines() if line.strip()}
        if len(hosts) > 1:
            raise ValueError('Use a single-node shared-memory PE for this workflow')
    if os.environ.get('PHENORADAR_TASK_MEM_MB'):
        memory_bytes = int(os.environ['PHENORADAR_TASK_MEM_MB']) * 1048576
    else:
        output = subprocess.check_output(['qstat', '-j', os.environ['JOB_ID']], text=True, timeout=10)
        match = re.search(r'(?:^|[,\s])(?:s_vmem|h_vmem|mem_req)=([0-9.]+)([kKmMgGtT]?)', output, re.M)
        if not match:
            raise ValueError('Cannot read finite AGE memory allocation; submit with an explicit memory request')
        power = {'': 0, 'k': 1, 'm': 2, 'g': 3, 't': 4}[match[2].lower()]
        memory_bytes = int(float(match[1]) * 1024 ** power) * slots
    budget = memory_bytes // 1000000 - 4000
    if slots < 1 or budget < 1:
        raise ValueError('Allocate positive CPUs and more than 4 GB for workflow overhead')
    return slots, budget


def resolve(settings):
    result = copy.deepcopy(settings)
    result['scheduler'] = detect(settings.get('scheduler', 'auto'))
    if result['scheduler'] == 'age':
        age = result.setdefault('age', {})
        if not age.get('parallel_environment') or not age.get('memory_resource'):
            pe, memory = age_capabilities()
            age['parallel_environment'] = age.get('parallel_environment') or pe
            age['memory_resource'] = age.get('memory_resource') or memory
        age.setdefault('extra_resources', {})
    return result


def parse_id(output, backend):
    text = output.strip()
    if backend == 'slurm':
        identifier = text.split(';')[0]
        if identifier.isdigit():
            return identifier
    else:
        match = re.fullmatch(r'(\d+)(?:\.\d+(?:-\d+(?::\d+)?)?)?', text)
        if match:
            return match[1]
    raise ValueError(f'Unrecognized {backend} submission response: {text!r}; reconcile the receipt before retrying')


def batch_command(settings, job, root, name, stdout, stderr, *, script=None,
                  indices=None, concurrency=None, afterok=(), afterany=()):
    """Keep resource units and dependencies consistent across both schedulers."""
    from phase_config import time_minutes
    backend = settings.get('scheduler', 'slurm')
    if backend == 'slurm':
        cmd = ['sbatch', '--parsable', '--nodes=1', '--ntasks=1', f"--cpus-per-task={job['cpus']}",
               f"--mem={job['mem_gb'] * 1000}M", f"--time={job['time']}", '--chdir=' + str(root),
               '--job-name=' + name, '--output=' + str(stdout), '--error=' + str(stderr)]
        if settings.get('partition'):
            cmd.append('--partition=' + settings['partition'])
        if indices:
            value = ','.join(map(str, indices))
            if concurrency is not None:
                value += '%' + str(concurrency)
            cmd.append('--array=' + value)
        dependencies = []
        if afterany:
            dependencies.append('afterany:' + ':'.join(map(str, afterany)))
        if afterok:
            dependencies.append('afterok:' + ':'.join(map(str, afterok)))
            cmd.append('--kill-on-invalid-dep=yes')
        if dependencies:
            cmd.append('--dependency=' + ','.join(dependencies))
    else:
        age = settings['age']
        if {age['memory_resource'], 'h_rt'} & age.get('extra_resources', {}).keys():
            raise ValueError('AGE extra_resources cannot override per-job memory or time')
        per_slot_mib = math.ceil(job['mem_gb'] * 1000 / job['cpus'])
        resources = {age['memory_resource']: f'{per_slot_mib}M',
                     'h_rt': str(time_minutes(job['time']) * 60), **age.get('extra_resources', {})}
        # % tokens in the common log template become AGE's runtime variables.
        def log_path(value):
            return str(value).replace('%A', '$JOB_ID').replace('%j', '$JOB_ID').replace('%a', '$TASK_ID')
        cmd = ['qsub', '-terse', '-cwd', '-V', '-S', '/bin/bash', '-N', name[:120],
               '-pe', age['parallel_environment'], str(job['cpus']), '-l',
               ','.join(str(k) if v is True else f'{k}={v}' for k, v in resources.items()),
               '-o', log_path(stdout), '-e', log_path(stderr), '-v',
               f"PHENORADAR_TASK_CPUS={job['cpus']},PHENORADAR_TASK_MEM_MB={per_slot_mib * job['cpus']}"]
        if age.get('queue'):
            cmd += ['-q', age['queue']]
        if indices:
            cmd += ['-t', f'1-{len(indices)}']
            if concurrency is not None:
                cmd += ['-tc', str(concurrency)]
        if afterok or afterany:
            cmd += ['-hold_jid', ','.join(dict.fromkeys(map(str, [*afterok, *afterany])))]
    return [*cmd, *([str(script)] if script is not None else [])]


def array_size(settings):
    if settings.get('scheduler', 'slurm') == 'slurm':
        from phase_config import resolve_array_size
        return resolve_array_size()
    output = subprocess.check_output(['qconf', '-sconf'], text=True, timeout=10)
    match = re.search(r'^max_aj_tasks\s+(\d+)\s*$', output, re.M)
    if not match or int(match[1]) < 1:
        raise ValueError('Cannot read AGE max_aj_tasks from qconf -sconf')
    return int(match[1])


def active_ids(backend):
    if backend == 'slurm':
        return set(subprocess.check_output(['squeue', '--noheader', '--user',
                   pwd.getpwuid(os.getuid()).pw_name, '--format=%i'], text=True).split())
    tree = ET.fromstring(subprocess.check_output(['qstat', '-xml'], text=True, timeout=10))
    return {item.findtext('JB_job_number') for item in tree.findall('.//job_list')}


def age_states(identifiers, run=None):
    """Never infer success from disappearance; require a matching qacct record."""
    run = run or subprocess.run
    result = {}
    query = run(['qstat', '-xml'], capture_output=True, text=True, timeout=10)
    if query.returncode:
        return result
    try:
        jobs = ET.fromstring(query.stdout).findall('.//job_list')
    except ET.ParseError:
        return result

    def matches(specification, task):
        if not specification or not task:
            return True
        for part in specification.split(','):
            match = re.fullmatch(r'(\d+)(?:-(\d+)(?::(\d+))?)?', part.strip())
            if match:
                first, last, step = int(match[1]), int(match[2] or match[1]), int(match[3] or 1)
                if first <= int(task) <= last and (int(task) - first) % step == 0:
                    return True
        return False

    for identifier in identifiers:
        base, _, task = identifier.partition('_')
        for node in jobs:
            if node.findtext('JB_job_number') == base and matches(node.findtext('tasks'), task):
                state = node.findtext('state', '')
                result[identifier] = {'state': 'FAILED' if 'E' in state else
                                      'RUNNING' if 'r' in state or 't' in state else 'PENDING',
                                      'exit_code': None, 'source': 'qstat'}
                break
    pending = set(identifiers) - result.keys()
    for base in sorted({identifier.split('_')[0] for identifier in pending}):
        response = run(['qacct', '-j', base], capture_output=True, text=True, timeout=10)
        if response.returncode:
            continue  # Accounting may lag; unknown is not success.
        for block in re.split(r'^=+\s*$', response.stdout, flags=re.M):
            fields = {}
            for line in block.splitlines():
                parts = line.split(None, 1)
                if len(parts) == 2:
                    fields[parts[0]] = parts[1].strip()
            if fields.get('jobnumber') != base or 'exit_status' not in fields or 'failed' not in fields:
                continue
            task = fields.get('taskid', '')
            identifier = base + ('_' + task if task.isdigit() and int(task) else '')
            if identifier in pending:
                success = fields['exit_status'] == '0' and fields['failed'].split()[0] == '0'
                result[identifier] = {'state': 'COMPLETED' if success else 'FAILED',
                                      'exit_code': '0:0' if success else fields['exit_status'] + ':1',
                                      'source': 'qacct'}
    return result


def storage_prefix(pipeline, root, usage, label):
    """All scheduler jobs use the same node-local storage and measurement wrapper."""
    return [sys.executable, str(Path(__file__).with_name('storage_runtime.py')), 'run',
            '--pipeline', str(pipeline), '--root', str(root), '--usage', str(usage), '--label', label, '--']


def storage_environment(pipeline, root, usage):
    return {'PHENORADAR_STORAGE_CONFIG': str(pipeline), 'PHENORADAR_PROJECT_ROOT': str(root),
            'PHENORADAR_USAGE_ROOT': str(usage), 'PHENORADAR_SUBMIT_PYTHON': sys.executable}


def worker_guard(root, owner):
    records = Path(root) / '.snakemake/cluster_submissions'
    pending = []
    for path in records.glob('*.json'):
        saved = json.loads(path.read_text())
        if saved.get('analysis') != str(owner):
            continue
        if saved['state'] != 'accepted':
            raise ValueError(f'Unresolved worker submission: {path}')
        pending.append(saved['job_id'])
    if pending and set(pending) & active_ids('age'):
        raise ValueError('analysis still has queued/running AGE workers')


def wrap_slurm(argv):
    """Wrap the Slurm plugin's --wrap command on the execution node, not here."""
    if not os.environ.get('PHENORADAR_STORAGE_CONFIG'):
        return argv
    prefix = storage_prefix(os.environ['PHENORADAR_STORAGE_CONFIG'], os.environ['PHENORADAR_PROJECT_ROOT'],
                            os.environ['PHENORADAR_USAGE_ROOT'], 'workflow')
    result = list(argv)
    for i, arg in enumerate(result):
        if arg == '--wrap':
            result[i + 1] = shlex.join([*prefix, 'bash', '-c', result[i + 1]])
            return result
        if arg.startswith('--wrap='):
            result[i] = '--wrap=' + shlex.join([*prefix, 'bash', '-c', arg.split('=', 1)[1]])
            return result
    raise ValueError('Unrecognized Slurm worker command; cannot safely apply common storage settings')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    command = sub.add_parser('detect')
    command.add_argument('--optional', action='store_true')
    command.add_argument('--scheduler', default='auto')
    sub.add_parser('allocation')
    command = sub.add_parser('profile')
    command.add_argument('--root', required=True)
    command.add_argument('--pipeline', required=True)
    command = sub.add_parser('submit-job')
    command.add_argument('--settings', required=True)
    command.add_argument('jobscript')
    command = sub.add_parser('status-job')
    command.add_argument('--state', required=True)
    command.add_argument('job_id')
    args = parser.parse_args()
    if args.command == 'detect':
        value = detect(args.scheduler, required=not args.optional)
        if value is None and all(shutil.which(c) for c in ('sbatch', 'squeue', 'sacct', 'qsub', 'qstat', 'qacct')):
            value = detect(args.scheduler)
        print(value or 'local')
    elif args.command == 'allocation':
        print(*allocation())
    elif args.command == 'profile':
        from phase_config import execution_config, read_yaml, validate_slurm, write_profile
        root = Path(args.root).resolve()
        execution = resolve(validate_slurm(execution_config(read_yaml(root / 'config/analysis.yaml')), build=False))
        directory = root / '.snakemake/profiles' / uuid.uuid4().hex
        print(write_profile(directory, execution, root=root, pipeline=args.pipeline,
                            usage=root / '.snakemake/storage'))
    elif args.command == 'status-job':
        status_file = Path(args.state) / (args.job_id + '.status.json')
        saved = json.loads(status_file.read_text())
        try:
            value = age_states([args.job_id]).get(args.job_id)
        except (OSError, subprocess.SubprocessError):
            value = None
        if value is None:
            missing = saved.get('missing_since') or time.time()
            if time.time() - missing > 600:
                raise ValueError('AGE accounting is unavailable for over 10 minutes; completion is unconfirmed')
            if not saved.get('missing_since'):
                write_json(status_file, dict(saved, missing_since=missing))
            print('running')
            return
        if saved.get('missing_since'):
            write_json(status_file, dict(saved, missing_since=None))
        print('success' if value['state'] == 'COMPLETED' else 'failed' if value['state'] == 'FAILED' else 'running')
    else:
        from snakemake.utils import read_job_properties
        from dataset_assets import locked, record
        config = json.loads(Path(args.settings).read_text())
        properties = read_job_properties(args.jobscript)
        resources = properties.get('resources', {})
        job = {'cpus': properties['threads'], 'mem_gb': max(1, math.ceil(resources.get('mem_mb', 8000) / 1000)),
               'time': str(int(resources.get('runtime', 1440)))}
        script = Path(args.jobscript).absolute()
        state = Path(config['root']) / '.snakemake/cluster_submissions'
        signature = hashlib.sha256((str(script) + script.read_text()).encode()).hexdigest()
        receipt = state / (signature + '.json')
        with locked(state / (signature + '.lock')):
            if receipt.exists():
                saved = json.loads(receipt.read_text())
                if saved['state'] != 'accepted':
                    raise ValueError(f'Unresolved AGE submission: {receipt}')
                print(saved['job_id'])
                return
            token = uuid.uuid4().hex
            wrapper = state / (signature + '.sh')
            prefix = storage_prefix(config['pipeline'], config['root'], config['usage'], properties.get('rule', 'group'))
            wrapper.write_text('#!/usr/bin/env bash\nset -euo pipefail\nexec ' + shlex.join([*prefix, 'bash', str(script)]) + '\n')
            wrapper.chmod(0o755)
            cmd = batch_command(config['execution'], job, config['root'], 'pr_' + token,
                                str(script) + '.%j.out', str(script) + '.%j.err', script=wrapper)
            saved = {'state': 'unknown', 'created_at': now(), 'command': cmd, 'script': record(script),
                     'analysis': os.environ.get('PHENORADAR_ANALYSIS_PATH')}
            write_json(receipt, saved)
            identifier = parse_id(subprocess.check_output(cmd, cwd=config['root'], text=True, timeout=60), 'age')
            saved.update(state='accepted', job_id=identifier)
            write_json(state / (identifier + '.status.json'), {'job_id': identifier, 'missing_since': None})
            write_json(receipt, saved)
            print(identifier)


if __name__ == '__main__':
    try:
        main()
    except (OSError, ValueError, subprocess.SubprocessError) as error:
        print(str(error), file=sys.stderr)
        raise SystemExit(1)
