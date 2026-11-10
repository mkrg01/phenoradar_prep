"""Run trimAl once and project its columns onto the original protein residues."""
from operator import itemgetter
from pathlib import Path
import re
import shutil
import subprocess
import string
import tempfile

from busco_phylogeny import fasta_records
from common import file_record

AMINO = set('ACDEFGHIKLMNPQRSTVWY')
MODES = {'gappyout', 'automated1'}
MISSING = str.maketrans({c: '-' for c in string.ascii_uppercase + '*.' if c not in AMINO})


def trim_records(records, command, mode, directory):
    if mode not in MODES:
        raise ValueError('trimAl mode must be gappyout or automated1')
    binary = shutil.which(command)
    if binary is None:
        raise FileNotFoundError(f'trimAl executable unavailable: {command}')
    binary = str(Path(binary).resolve())
    names = [name for name, _ in records]
    widths = {len(seq) for _, seq in records}
    if not records or len(set(names)) != len(names) or len(widths) != 1 or 0 in widths:
        raise ValueError('alignment must have unique names and equal nonempty sequence lengths')
    viewed = {name: seq.translate(MISSING) for name, seq in records}
    Path(directory).mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.trimal-', dir=directory) as temporary:
        source, output = Path(temporary) / 'missing_as_gaps.faa', Path(temporary) / 'trimmed.faa'
        with source.open('w') as handle:
            for name in names:
                handle.write(f'>{name}\n{viewed[name]}\n')
        argv = [binary, '-in', str(source), '-out', str(output), '-fasta', '-keepheader',
                '-keepseqs', '-colnumbering', '-' + mode]
        result = subprocess.run(argv, check=True, stdout=subprocess.PIPE, text=True)
        maps = re.findall(r'^#ColumnsMap\t(.*)$', result.stdout, re.MULTILINE)
        if len(maps) != 1 or not re.fullmatch(r'\d+(?:,\s*\d+)*\s*', maps[0]):
            raise ValueError('trimAl did not return a valid column map')
        selected = [int(value) for value in maps[0].split(',')]
        if selected != sorted(set(selected)) or selected[-1] >= next(iter(widths)):
            raise ValueError('trimAl column map is duplicated, unordered or out of range')
        project = itemgetter(*selected)
        trimmed = list(fasta_records(output))
        if len(trimmed) != len(names) or {n for n, _ in trimmed} != set(names):
            raise ValueError('trimAl changed the species set despite -keepseqs')
        if any(seq != ''.join(project(viewed[name])) for name, seq in trimmed):
            raise ValueError('trimAl residues disagree with its column map')
    projected = [(name, ''.join(project(seq))) for name, seq in records]
    return projected, selected, argv, file_record(binary)
