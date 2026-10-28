"""Small, per-analysis execution observations; never a completion certificate."""
import json
import logging
import os
from pathlib import Path
import re
import time
import uuid

from common import now, write_json


def read_state(path):
    try:
        path = Path(path)
        if path.stat().st_size > 65536:
            return {}
        value = json.loads(path.read_text())
        return value if isinstance(value, dict) else {}
    except (OSError, ValueError):
        return {}


class RunState:
    def __init__(self, analysis, target, execution=None):
        self.path = Path(analysis) / 'jobs/run_state.json'
        self.token = uuid.uuid4().hex
        self.data = {'schema_version': 1, 'run_token': self.token, 'target': target,
                     'controller_id': os.environ.get('SLURM_JOB_ID'),
                     'execution': str(Path(execution).resolve()) if execution else None,
                     'started_at': now(), 'phase': 'starting'}

    def phase(self, phase, **fields):
        current = read_state(self.path)
        if current.get('run_token') == self.token:
            self.data.update(current)
        if self.data['phase'] != phase:
            self.data.pop('progress', None)
        self.data.update(fields, phase=phase, updated_at=now())
        try:
            write_json(self.path, self.data)
        except OSError:
            pass  # Execution/completion receipts remain authoritative if observation storage fails.


class ProgressHandler(logging.Handler):
    def __init__(self, path, token):
        super().__init__()
        self.path, self.token, self.last_write = Path(path), token, 0

    def emit(self, record):
        try:
            event = getattr(record, 'event', None)
            event = getattr(event, 'value', event)
            submitted = re.search(r'has been submitted with SLURM jobid (\d+)', record.getMessage())
            progress = event == 'progress' and hasattr(record, 'done') and hasattr(record, 'total')
            if not (progress or submitted or record.levelno >= logging.ERROR):
                return
            state = read_state(self.path)
            if state.get('run_token') != self.token:
                return
            if state.get('controller_id') != os.environ.get('SLURM_JOB_ID'):
                return  # Worker Snakefiles must not overwrite controller observations.
            if progress:
                if time.monotonic() - self.last_write < 5 and record.done < record.total:
                    return
                state['progress'] = {'done': record.done, 'total': record.total,
                                     'unit': 'workflow_steps'}
            if submitted:
                ids = state.get('recent_worker_ids', [])
                state['recent_worker_ids'] = list(dict.fromkeys([*ids, submitted[1]]))[-64:]
            if record.levelno >= logging.ERROR:
                state['last_error'] = record.getMessage()[-1600:]
            state['updated_at'] = now()
            write_json(self.path, state)
            self.last_write = time.monotonic()
        except Exception:
            pass  # Observation failures must not change scientific execution.


def install_progress(workflow):
    token, path = os.environ.get('PHENORADAR_RUN_TOKEN'), os.environ.get('PHENORADAR_RUN_STATE')
    if not token or not path or os.environ.get('PHENORADAR_READ_ONLY'):
        return
    from snakemake_interface_executor_plugins.settings import ExecMode
    if workflow.exec_mode != ExecMode.DEFAULT:
        return
    from snakemake.logging import logger
    logger.addHandler(ProgressHandler(path, token))
