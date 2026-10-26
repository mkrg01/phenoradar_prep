"""Rate-limited sbatch with reconciliation of ambiguous timeout responses."""
from contextlib import contextmanager
import datetime
import fcntl
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import uuid


class UnknownSubmission(RuntimeError):
    pass


def check_analysis_submissions(state_dir, analysis, queued):
    """A restarted controller must not overlap surviving or ambiguous workers."""
    for path in Path(state_dir).glob('*.json'):
        state = json.loads(path.read_text())
        if state.get('analysis') != str(Path(analysis).resolve()):
            continue
        if state['state'] in {'submitting', 'unknown'}:
            raise ValueError(f'unresolved worker submission in {path}; reconcile before restarting')
        if state.get('job_id') in queued:
            raise ValueError(f'analysis still has queued/running worker {state["job_id"]}')


def option(argv, name, default=''):
    for i, value in enumerate(argv):
        if value.startswith(name + '='):
            return value.split('=', 1)[1]
        if value == name and i + 1 < len(argv):
            return argv[i + 1]
    return default


def save(path, data):
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(data, sort_keys=True) + '\n')
    temporary.replace(path)


@contextmanager
def locked(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('a') as handle:
        fcntl.flock(handle, fcntl.LOCK_EX)
        yield


def reconcile(token, started, workdir, run=subprocess.run):
    """Positive identification only. An empty/unavailable queue is not rejection."""
    comment = 'phenoradar:' + token
    start = datetime.datetime.fromtimestamp(started - 120).strftime('%Y-%m-%dT%H:%M:%S')
    commands = [
        ['squeue', '--noheader', '--user', str(os.getuid()), '--format=%i|%k|%Z'],
        ['sacct', '-X', '--noheader', '--parsable2', '--user', str(os.getuid()),
         '--starttime', start, '--format=JobIDRaw,Comment%256,WorkDir%1024'],
    ]
    found = set()
    for command in commands:
        try:
            result = run(command, capture_output=True, text=True, timeout=20)
        except (OSError, subprocess.TimeoutExpired):
            continue
        if result.returncode:
            continue
        for line in result.stdout.splitlines():
            fields = line.strip().split('|')
            if len(fields) >= 3 and comment in fields[1].strip().split(';') and fields[0].isdigit():
                if Path(fields[2].strip()).resolve() == Path(workdir).resolve():
                    found.add(fields[0])
    if len(found) > 1:
        raise UnknownSubmission(f'multiple Slurm jobs match submission token {token}')
    return next(iter(found), None)


def submit(argv, real_sbatch, state_dir, run=subprocess.run, sleep=time.sleep, clock=time.time, stdin=None):
    state_dir = Path(state_dir)
    state_dir.mkdir(parents=True, exist_ok=True)
    workdir = option(argv, '--chdir', option(argv, '-D', os.getcwd()))
    identity = [option(argv, '--job-name'), option(argv, '--comment'), workdir]
    key = hashlib.sha256(json.dumps(identity).encode()).hexdigest()
    signature = hashlib.sha256(json.dumps([argv, stdin]).encode()).hexdigest()
    receipt = state_dir / f'{key}.json'
    with locked(state_dir / f'{key}.lock'):
        previous = json.loads(receipt.read_text()) if receipt.exists() else None
        if previous and previous['state'] == 'accepted' and previous['signature'] == signature:
            return previous['job_id']
        if previous and previous['state'] in {'submitting', 'unknown'}:
            state = previous
        else:
            state = {'token': uuid.uuid4().hex, 'signature': signature,
                     'started': clock(), 'state': 'new', 'attempts': 0,
                     'analysis': os.environ.get('PHENORADAR_ANALYSIS_PATH')}

        def recover():
            for delay in (5, 10, 20):
                sleep(delay)
                identifier = reconcile(state['token'], state['started'], workdir, run)
                if identifier:
                    state.update(state='accepted', job_id=identifier)
                    save(receipt, state)
                    return identifier
            state['state'] = 'unknown'
            save(receipt, state)
            raise UnknownSubmission(f'Unresolved Slurm submission; no automatic resubmission. Inspect {receipt}')

        if state['state'] in {'submitting', 'unknown'}:
            return recover()
        original_comment = option(argv, '--comment')
        comment = (original_comment + ';' if original_comment else '') + 'phenoradar:' + state['token']
        clean, i = [], 0
        while i < len(argv):
            value = argv[i]
            if value == '--comment':
                i += 2
                continue
            if not value.startswith('--comment='):
                clean.append(value)
            i += 1
        call = [real_sbatch, '--comment=' + comment, *clean]
        for attempt in range(3):
            with locked(state_dir / 'rate.lock'):
                gate = state_dir / 'rate.json'
                last = json.loads(gate.read_text())['time'] if gate.exists() else 0
                sleep(max(0, 2 - (clock() - last)))
                save(gate, {'time': clock()})
            state.update(state='submitting', attempts=attempt + 1)
            save(receipt, state)
            try:
                result = run(call, capture_output=True, text=True, **({'input': stdin} if stdin is not None else {}))
            except OSError:
                state['state'] = 'rejected'; save(receipt, state)
                raise
            identifier = result.stdout.strip().split(';')[0]
            if result.returncode == 0 and identifier.isdigit():
                state.update(state='accepted', job_id=identifier)
                save(receipt, state)
                return identifier
            message = result.stderr.strip()
            state['error'] = message
            if 'Slurm temporarily unable to accept job' in message:
                state['state'] = 'rejected'; save(receipt, state)
                if attempt < 2:
                    sleep(5 * 2 ** attempt)
                    continue
                raise RuntimeError(message)
            rejected = ('invalid partition', 'invalid account', 'invalid qos',
                        'invalid option', 'unrecognized option',
                        'requested node configuration is not available',
                        'job violates accounting/qos policy', 'access/permission denied')
            if result.returncode != 0 and any(x in message.lower() for x in rejected):
                state['state'] = 'rejected'; save(receipt, state)
                raise RuntimeError(message)
            return recover()


if __name__ == '__main__':
    try:
        argv = sys.argv[1:]
        script = sys.stdin.read() if any(x in {'-', '/dev/stdin'} for x in argv) else None
        print(submit(argv, os.environ['PHENORADAR_REAL_SBATCH'],
                     os.environ['PHENORADAR_SUBMIT_STATE'], stdin=script))
    except (OSError, RuntimeError, KeyError, ValueError) as error:
        print(f'sbatch: {error}', file=sys.stderr)
        raise SystemExit(75 if isinstance(error, UnknownSubmission) else 1)
