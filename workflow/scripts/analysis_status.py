"""Read-only, bounded observations of analysis execution and completion."""
import datetime
import json
import os
from pathlib import Path
import re
import subprocess

from analysis_runtime import read_state
from common import now
from layout import SEQUENCE_EXPRESSION, CONSUMER_METADATA, BASE_METADATA

FAILED = {'FAILED', 'OUT_OF_MEMORY', 'TIMEOUT', 'CANCELLED', 'NODE_FAIL',
          'BOOT_FAIL', 'DEADLINE', 'PREEMPTED', 'REVOKED'}
ACTIVE = {'RUNNING', 'COMPLETING', 'CONFIGURING', 'SUSPENDED', 'RESIZING', 'STAGE_OUT'}


def submissions(path):
    return [(p, json.loads(p.read_text())) for p in sorted((Path(path) / 'jobs').glob('submission_*.json'))
            if re.fullmatch(r'submission_[0-9]+\.json', p.name)]


def timestamp(value):
    try:
        return datetime.datetime.fromisoformat(value.replace('Z', '+00:00')).timestamp()
    except (ValueError, TypeError, AttributeError):
        return None


def header_time(path):
    try:
        with Path(path).open() as handle:
            prefix = handle.read(8192)
        match = re.search(r'"created_at"\s*:\s*"([^"]+)"', prefix)
        return match[1] if match else None
    except OSError:
        return None


def observe_slurm(identifier, root, created_at, workers=(), run=None):
    run = run or subprocess.run
    result = {'job_id': identifier, 'state': 'unknown', 'reason': None,
              'exit_code': None, 'source': None, 'workers': [], 'errors': []}
    if not identifier or not str(identifier).isdigit():
        return result
    ids = {str(identifier), *(str(i) for i in workers if str(i).isdigit())}

    def query(command):
        try:
            response = run(command, capture_output=True, text=True, timeout=8)
            if response.returncode:
                result['errors'].append(response.stderr.strip()[:500] or 'Slurm query failed')
                return ''
            return response.stdout
        except (OSError, subprocess.TimeoutExpired) as error:
            result['errors'].append(str(error)[:500])
            return ''

    queue = query(['squeue', '--noheader', '--user', str(os.getuid()), '--format=%i|%T|%R|%Z'])
    for line in queue.splitlines():
        fields = line.strip().split('|')
        if len(fields) != 4 or fields[0] not in ids or fields[3] != str(Path(root).resolve()):
            continue
        item = {'job_id': fields[0], 'state': fields[1], 'reason': fields[2], 'source': 'squeue'}
        if fields[0] == str(identifier):
            result.update(item)
        else:
            result['workers'].append(item)
    if result['source']:
        return result
    start = timestamp(created_at)
    if start is None:
        result['errors'].append('Submission date unavailable; accounting was not queried')
        return result
    start = datetime.datetime.fromtimestamp(start - 86400).strftime('%Y-%m-%dT%H:%M:%S')
    history = query(['sacct', '-X', '--noheader', '--parsable2', '--jobs', str(identifier),
                     '--starttime', start, '--format=JobIDRaw,State,ExitCode,WorkDir%1024'])
    for line in history.splitlines():
        fields = line.strip().split('|')
        if len(fields) >= 4 and fields[0] == str(identifier) and fields[3] == str(Path(root).resolve()):
            result.update(state=fields[1].split()[0].rstrip('+'), exit_code=fields[2] or None, source='sacct')
            break
    return result


def log_paths(path, receipt):
    number = receipt[0].stem.removeprefix('submission_')
    identifier = receipt[1].get('job_id')
    return [Path(path) / f'jobs/logs/{number}_{identifier}.{suffix}' for suffix in ('out', 'err')]


def legacy_progress(path, receipt):
    """Old runs have no observation record. Read at most 64 KiB, never whole logs."""
    logs = log_paths(path, receipt)
    result = {'phase': None, 'progress': None, 'recent_worker_ids': []}
    for log in logs:
        try:
            with log.open('rb') as handle:
                handle.seek(max(0, log.stat().st_size - 65536))
                tail = handle.read(65536).decode(errors='replace')
        except OSError:
            continue
        found = list(re.finditer(r'(\d+) of (\d+) steps \(', tail))
        if found:
            result['progress'] = dict(zip(('done', 'total'), map(int, found[-1].groups())), unit='workflow_steps')
        result['recent_worker_ids'].extend(re.findall(r'has been submitted with SLURM jobid (\d+)', tail)[-64:])
        # Retain observations of old launchers without requiring their bundle.
        if 'Rule: phenoradar_inputs' in tail or 'rule phenoradar_inputs:' in tail:
            result['phase'] = 'exporting'
    return result


def required_outputs(pipeline, implementation=()):
    """Bounded presence checks for final branch markers, not per-OG validation."""
    required = ['run.json', 'metadata/samples.tsv', 'metadata/selection.json',
                'metadata/species_metadata.tsv', 'orthogroups/expression/tpm.tsv',
                'orthogroups/expression/mapping_qc.tsv']
    # Frozen analyses predating this producer retain their original contract.
    if any(Path(entry['path']).name == 'sequence_expression.py' for entry in implementation):
        required.append(SEQUENCE_EXPRESSION)
    if any(Path(entry['path']).name == 'consumer_metadata.py' for entry in implementation):
        prefix = 'filtered/' if pipeline.get('exclude_species') else ''
        required += [prefix + name for name in [BASE_METADATA, CONSUMER_METADATA,
                     str(Path(CONSUMER_METADATA).with_suffix('.json'))] if prefix + name not in required]
    if pipeline.get('alignment', {}).get('enabled'):
        required += ['orthogroups/alignments/provenance.json']
    if pipeline.get('kegg', {}).get('enabled'):
        required += ['kegg/provenance.json', 'kegg/annotation_provenance.json',
                     'kegg/ko_tpm_sum.tsv']
    phy = pipeline.get('phylogeny', {})
    for branch in phy.get('trees', []):
        prefix = f'phylogeny/{branch}'
        required += [f'{prefix}/species_tree.nwk', f'{prefix}/species_tree.json']
        for option, marker in [('dating', 'dating/provenance.json'),
                               ('taxonomy_check', 'taxonomy_check/summary.json'),
                               ('contrast_pairs', 'contrast/summary.json')]:
            if phy.get(option, {}).get('enabled'):
                required.append(f'{prefix}/{marker}')
    return required


def status(path, manifest=None, run=None):
    path = Path(path).resolve()
    manifest = manifest or json.loads((path / 'analysis.json').read_text())
    receipts = submissions(path)
    latest = receipts[-1] if receipts else None
    saved = latest[1] if latest else {}
    runtime = read_state(path / 'jobs/run_state.json')
    identifier = saved.get('job_id')
    runtime_matches = bool(runtime) and (
        str(runtime.get('controller_id') or '') == str(identifier or '') and
        runtime.get('target') == saved.get('target', runtime.get('target')) and
        (not latest or runtime.get('execution') == str(latest[0].with_suffix('.resources.json'))))
    if not runtime_matches:
        runtime = {}
    legacy = legacy_progress(path, latest) if latest and not runtime else {}
    phase = runtime.get('phase') or legacy.get('phase')
    controller = observe_slurm(identifier, manifest['root'], saved.get('created_at'),
                               runtime.get('recent_worker_ids', legacy.get('recent_worker_ids', [])), run=run)
    completion = path / 'completed.json'
    created_at = header_time(completion)
    started = timestamp(saved.get('created_at'))
    finished = timestamp(created_at)
    current = completion.is_file() and finished is not None and (not latest or started is not None and finished >= started)
    required = required_outputs(manifest.get('pipeline', {}), manifest.get('implementation', []))
    missing = [name for name in required if not (path / name).is_file()] if current else []
    complete = current and not missing and saved.get('target', 'all') == 'all'
    queue_state = controller['state']
    if saved.get('state') in {'submitting', 'unknown'}:
        state, action = 'unknown', 'Reconcile the unresolved submission before retrying.'
    elif saved.get('state') == 'rejected' or queue_state in FAILED or controller.get('exit_code') not in {None, '0:0'}:
        state, action = 'failed', 'Review the controller log and failure reason before submitting again.'
    elif queue_state == 'PENDING' and controller.get('reason') == 'DependencyNeverSatisfied':
        state, action = 'blocked', 'The dependency cannot succeed; inspect prerequisite jobs before retrying.'
    elif queue_state == 'PENDING':
        state, action = 'pending', 'Wait for the reported Slurm scheduling reason to clear.'
    elif queue_state in ACTIVE:
        state = 'finalizing' if phase in {'validating', 'exporting', 'recording_completion', 'complete'} else 'running'
        action = 'Wait for the current phase to finish.'
    elif complete and (queue_state == 'COMPLETED' or runtime.get('phase') == 'complete'):
        state, action = 'complete', 'Analysis outputs are ready in their native directories.'
    elif runtime.get('phase') in {'failed', 'interrupted'}:
        state, action = 'failed', 'Review the recorded error before submitting again.'
    elif current and missing:
        state, action = 'incomplete', 'Required outputs are missing; inspect submit --dry-run before recovery.'
    elif runtime.get('phase') == 'target_complete':
        state, action = 'target_complete', 'The requested target finished; whole-analysis completion is separate.'
    elif queue_state == 'COMPLETED':
        state = 'target_complete' if saved.get('target') != 'all' else 'incomplete'
        action = 'The requested target finished.' if state == 'target_complete' else 'Controller exited, but whole-analysis completion is not confirmed.'
    elif latest or runtime or completion.exists():
        state, action = 'unknown', 'Execution cannot be confirmed from the available observations.'
    else:
        state, action = 'prepared', 'Use submit --dry-run to inspect the pending work.'
    return {'schema_version': 2, 'name': manifest['name'], 'build': manifest['build']['path'],
            'state': state, 'phase': phase, 'observed_at': now(), 'next_action': action,
            'selection': manifest.get('selection', []), 'submissions': [r for _, r in receipts],
            'controller': controller, 'progress': runtime.get('progress') or legacy.get('progress'),
            'runtime': runtime, 'logs': [str(p) for p in log_paths(path, latest)] if latest else [],
            'completion': {'record': str(completion), 'created_at': created_at, 'current': bool(current),
                           'required_outputs_present': bool(complete), 'missing': missing},
            'verification': {'performed': False, 'scope': 'execution records and required output presence'},
            'storage': {'state': 'not_checked', 'counts': {}, 'jobs': [], 'apply': False}}


def print_status(report):
    controller = report['controller']
    print(f"Analysis: {report['name']}")
    print(f"State: {report['state']}" + (f" ({report['phase']})" if report['phase'] else ''))
    print(f"Controller: {controller['job_id'] or '-'} {controller['state']}")
    if controller.get('reason'):
        print('Slurm reason: ' + controller['reason'])
    progress = report.get('progress')
    if progress:
        print(f"Phase progress: {progress['done']}/{progress['total']} workflow steps")
    print('Verification: ' + ('passed' if report['verification']['performed'] else 'not requested (use --verify)'))
    print('Next: ' + report['next_action'])
    if report['runtime'].get('last_error'):
        print('Last recorded error: ' + report['runtime']['last_error'])
    for path in report['logs']:
        print('Log: ' + path)
    for error in controller.get('errors', []):
        print('Observation: ' + error)
