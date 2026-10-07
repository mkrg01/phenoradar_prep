"""Read GeneGalleon's own attempt evidence; reuse decisions stay in GeneGalleon."""
import json
import re
import subprocess
import sys
from pathlib import Path


STAGE_STEPS = {
    'assembly': {'transcriptome_longest_cds'},
    'busco': {'transcriptome_busco_longest_cds'},
    'quant': {'transcriptome_quant', 'transcriptome_merge'},
}


def read_json(path):
    if path.is_symlink() or path.stat().st_size > 8 * 1024 * 1024:
        raise ValueError('invalid GeneGalleon observation')
    return json.loads(path.read_text())


def observation(work):
    """Report recorded progress cheaply, without claiming fresh artifact verification."""
    directory = Path(work) / 'output/observations'
    attempts = []
    if directory.is_dir() and not directory.is_symlink():
        for attempt in directory.iterdir():
            try:
                run = read_json(attempt / 'run.json')
                if (isinstance(run, dict) and not attempt.is_symlink()
                        and run.get('schema') == 'genegalleon-observation-v1'
                        and run.get('attempt_id') == attempt.name
                        and run.get('workflow') == 'gg_transcriptome_generation'
                        and type(run.get('started_at_ns')) is int):
                    attempts.append((run['started_at_ns'], attempt, run))
            except (OSError, ValueError, KeyError, TypeError):
                continue
    if not attempts:
        return {'state': 'unavailable', 'completed_stages': [], 'current_stage': None}
    _, attempt, run = max(attempts, key=lambda row: row[0])
    receipts = []
    for path in attempt.glob('contract-*.json'):
        try:
            receipt = read_json(path)
            if (isinstance(receipt, dict) and receipt.get('schema') == 'genegalleon-contract-observation-v1'
                    and receipt.get('attempt_id') == attempt.name
                    and isinstance(receipt.get('step'), str)
                    and type(receipt.get('observed_at_ns')) is int):
                receipts.append(receipt)
        except (OSError, ValueError, TypeError):
            continue
    done = {r['step'] for r in receipts if r.get('manifest_sha256') and
            ((r.get('operation') == 'record' and r.get('exit_code') == 0) or
             (r.get('operation') == 'needs-run' and r.get('exit_code') == 1))}
    latest = max(receipts, key=lambda r: r.get('observed_at_ns', 0), default={})
    step = latest.get('step')
    current = next((stage for stage, steps in STAGE_STEPS.items() if step in steps), None)
    if step and current is None:
        current = 'assembly' if step.startswith('transcriptome_') else None
    return {'state': run.get('execution_state', 'unknown'), 'attempt': str(attempt),
            'completed_stages': [stage for stage, steps in STAGE_STEPS.items() if steps <= done],
            'current_stage': current, 'step': step, 'exit_code': run.get('exit_code'),
            'coverage': 'recorded-attempt-only'}


def verified_stages(work, repository):
    """Ask the pinned native API before publishing checkpoints from a failed run."""
    observed = observation(work)
    api = Path(repository) / 'workflow/support/workflow_api.py'
    if not observed.get('attempt') or not api.is_file():
        return set(), observed
    try:
        result = subprocess.run([sys.executable, '-B', str(api), 'preflight', '--attempt',
                                 observed['attempt'], '--workspace-root', str(work)],
                                capture_output=True, text=True, timeout=60, check=True)
        report = json.loads(result.stdout)
        if (not isinstance(report, dict) or report.get('schema') != 'genegalleon-api-v1'
                or report.get('command') != 'preflight'):
            raise ValueError('unsupported GeneGalleon API response')
        current = {r['step'] for r in report['contracts'] if r['state'] == 'verified_current'}
        observed['contracts'] = report['contracts']
        return {stage for stage, steps in STAGE_STEPS.items() if steps <= current}, observed
    except (OSError, ValueError, KeyError, TypeError, subprocess.SubprocessError) as error:
        observed['verification_error'] = str(error)
        return set(), observed


def scheduler_states(jobs):
    """Accounting evidence supplements native records after OOM/SIGKILL/timeouts."""
    ids = sorted({str(job.get('slurm_id') or job.get('job_id')) for job in jobs
                  if re.fullmatch(r'[0-9]+(?:_[0-9]+)?', str(job.get('slurm_id') or job.get('job_id')))})
    states = {}
    for offset in range(0, len(ids), 512):
        batch = ids[offset:offset + 512]
        try:
            result = subprocess.run(['sacct', '--noheader', '--parsable2', '--allocations',
                                     '--jobs', ','.join(batch), '--format=JobID%64,State%32,ExitCode'],
                                    capture_output=True, text=True, timeout=5, check=True)
            for line in result.stdout.splitlines():
                fields = line.split('|')
                state = fields[1].strip().split() if len(fields) >= 3 else []
                if state and fields[0].strip() in batch:
                    states[fields[0].strip()] = {'state': state[0],
                                                 'exit_code': fields[2].strip()}
        except (OSError, subprocess.SubprocessError):
            continue
    return states
