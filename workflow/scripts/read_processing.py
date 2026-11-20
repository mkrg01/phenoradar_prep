"""Summarize read acquisition and fastp evidence before native diagnostics retire."""
import csv
import json
import math
from pathlib import Path


FIELDS = ('read_acquisition_method', 'fastp_applied')
METHODS = {'amalgkit': 'amalgkit', 'public-original': 'genegalleon_fallback'}


def _number(value):
    try:
        number = float(value)
        return number if math.isfinite(number) and number >= 0 and number.is_integer() else None
    except (ValueError, TypeError):
        return None


def _fastp_stats(path, run):
    if not path.is_file() or path.is_symlink() or path.stat().st_size > 8 * 1024 * 1024:
        return 'unknown'
    with path.open(encoding='utf-8', newline='') as handle:
        rows = [row for row in csv.DictReader(handle, delimiter='\t') if row.get('run') == run]
    if len(rows) != 1:
        return 'unknown'
    count = _number(rows[0].get('num_fastp_in'))
    if count is not None and count > 0:
        return 'yes'
    written = _number(rows[0].get('num_written'))
    if count == 0 and written is not None and written > 0:
        return 'no'
    return 'unknown'


def collect(output, native, run):
    """Use completion records, not requested settings or fallback bp_fastp placeholders."""
    result = dict.fromkeys(FIELDS, 'unknown')
    root = Path(output) / 'amalgkit_getfastq' / native
    path = root / 'getfastq_completion.json'
    try:
        if not path.is_file() or path.is_symlink() or path.stat().st_size > 8 * 1024 * 1024:
            return result
        data = json.loads(path.read_text(encoding='utf-8'))
        if not isinstance(data, dict) or data.get('status') != 'complete':
            return result
        entries = data.get('runs')
        if not isinstance(entries, list):
            return result
        matches = [entry for entry in entries if isinstance(entry, dict) and entry.get('run') == run]
        if len(matches) != 1 or matches[0].get('status', 'complete') != 'complete':
            return result
        entry = matches[0]
        source = data.get('read_source')
        result['read_acquisition_method'] = METHODS.get(source, 'unknown')
        # A future producer can report preprocessing independently of acquisition.
        if 'fastp_applied' in entry or 'fastp_applied' in data:
            applied = entry.get('fastp_applied', data.get('fastp_applied'))
            if type(applied) is bool:
                result['fastp_applied'] = 'yes' if applied else 'no'
            elif isinstance(applied, str) and applied in {'yes', 'no', 'unknown'}:
                result['fastp_applied'] = applied
        else:
            result['fastp_applied'] = _fastp_stats(root / run / 'getfastq_stats.tsv', run)
            if (result['fastp_applied'] == 'unknown' and source == 'public-original'
                    and data.get('schema_version') in {3, 4}):
                # These native schemas describe raw rescue FASTQs. Their bp_fastp
                # fields repeat input bases even though fastp was not executed.
                result['fastp_applied'] = 'no'
    except (OSError, ValueError, TypeError, csv.Error):
        # Missing or unsupported historical evidence must not invent a method.
        pass
    return result


def summary(processing):
    """One row covers assembly and quantification; keep stage differences explicit."""
    stages = processing.get('read_processing', {})
    result = {}
    for field in FIELDS:
        values = [stages.get(stage, {}).get(field, 'unknown') for stage in ('assembly', 'quant')]
        result[field] = ('unknown' if 'unknown' in values else
                         values[0] if values[0] == values[1] else 'mixed')
    return result
