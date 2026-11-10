import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

from align_orthogroups import align, fasta_records
from alignment_batches import plan, run_batch, batch_targets
from alignment_trimming import trim_records
from configuration import validate_keys


def trimmer(tmp_path, columns=(0, 2, 3)):
    script = tmp_path / 'trimal'
    script.write_text(f'''#!{sys.executable}
import sys
from pathlib import Path
args = sys.argv
records = []
for line in Path(args[args.index('-in') + 1]).read_text().splitlines():
    if line.startswith('>'):
        records.append([line, ''])
    else:
        records[-1][1] += line
columns = {list(columns)!r}
Path(args[args.index('-out') + 1]).write_text(''.join(name + '\\n' + ''.join(seq[i] for i in columns) + '\\n' for name, seq in records))
print('#ColumnsMap\\t' + ', '.join(map(str, columns)))
''')
    script.chmod(0o755)
    return str(script)


def test_projection_restores_unknowns_and_retains_all_gap_copies(tmp_path):
    records = [('Plant_g1', 'MAX*'), ('Plant_g2', '----')]
    projected, columns, _, _ = trim_records(records, trimmer(tmp_path), 'gappyout', tmp_path)
    assert projected == [('Plant_g1', 'MX*'), ('Plant_g2', '---')]
    assert columns == [0, 2, 3]


@pytest.mark.parametrize('columns', [(0, 0), (2, 1), (0, 9)])
def test_invalid_column_maps_fail_before_publication(tmp_path, columns):
    # A deliberately broken trimmer reports invalid maps without indexing them.
    command = tmp_path / 'bad'
    command.write_text(f'#!{sys.executable}\nprint("#ColumnsMap\\t' + ', '.join(map(str, columns)) + '")\n')
    command.chmod(0o755)
    with pytest.raises(ValueError, match='column map'):
        trim_records([('Plant_g1', 'MAAA')], str(command), 'gappyout', tmp_path)


def test_batch_trims_legacy_once_and_mode_change_invalidates_receipts(tmp_path):
    inputs = tmp_path / 'inputs'; inputs.mkdir()
    fasta = inputs / 'OG1.faa'; fasta.write_text('>Plant_g1\nMACD\n')
    (inputs / 'provenance.json').write_text(json.dumps({'orthogroups': ['OG1']}))
    out, reports, plans = [tmp_path / name for name in ['out', 'reports', 'plans']]
    align(fasta, out / fasta.name, reports / 'OG1.json', command='unavailable-famsa')
    original = json.loads((reports / 'OG1.json').read_text())['input']
    command = trimmer(tmp_path)
    receipt = out / 'batches/batch_000000.json'
    for mode in ['gappyout', 'automated1']:
        plan(inputs, plans, trimal_mode=mode)
        assert not Path(batch_targets(plans, out, reports)[0]).exists()
        run_batch(plans / 'batch_000000.json', inputs, out, reports, receipt,
                  command='unavailable-famsa', trimal_mode=mode, trimal_command=command)
        assert dict(fasta_records(out / fasta.name)) == {'Plant_g1': 'MCD'}
        report = json.loads((reports / 'OG1.json').read_text())
        assert report['input'] == original
        assert report['trimming']['mode'] == mode
        assert report['trimming']['raw_sites'] == 4
        saved = (out / fasta.name).stat().st_mtime_ns
        run_batch(plans / 'batch_000000.json', inputs, out, reports, receipt,
                  command='unavailable-famsa', trimal_mode=mode, trimal_command=command)
        assert json.loads(receipt.read_text())['reused'] == 1
        assert (out / fasta.name).stat().st_mtime_ns == saved


def test_failed_trimming_preserves_existing_alignment(tmp_path):
    fasta = tmp_path / 'input.faa'; fasta.write_text('>Plant_g1\nMACD\n')
    output = tmp_path / 'out.faa'; output.write_text('old result\n')
    with pytest.raises(FileNotFoundError, match='trimAl'):
        align(fasta, output, tmp_path / 'report.json', trimal_mode='gappyout',
              trimal_command='unavailable-trimal')
    assert output.read_text() == 'old result\n'
    assert not (tmp_path / 'report.json').exists()


def test_og_trim_mode_is_validated_before_scheduling():
    with pytest.raises(ValueError, match='alignment.trimal_mode'):
        validate_keys({'alignment': {'trimal_mode': 'typo'}})


def test_real_batch_workflow_runs_alignment_and_trim_in_one_job(workflow_project):
    snakemake = os.environ.get('SNAKEMAKE_BIN') or shutil.which('snakemake')
    trimal = os.environ.get('TRIMAL_BIN') or shutil.which('trimal')
    if not snakemake or not trimal or not shutil.which('famsa'):
        pytest.skip('Snakemake, FAMSA and trimAl required')
    root = Path(__file__).resolve().parents[1]
    inputs = workflow_project / 'work/test/orthogroups/alignments/selected'
    inputs.mkdir(parents=True)
    (inputs / 'OG1.faa').write_text('>Plant_g1\nMACD\n>Plant_g2\nMAD\n')
    (inputs / 'provenance.json').write_text(json.dumps({'orthogroups': ['OG1']}))
    plan(inputs, inputs.parent / 'batches', trimal_mode='gappyout')
    # Exercise the production worker rule after its planning checkpoint has
    # completed, independently of upstream ODB/database preparation.
    rules = (root / 'workflow/rules/alignments.smk').read_text()
    worker = 'rule align_orthogroup_batch:' + rules.split('rule align_orthogroup_batch:', 1)[1].split('\n\nrule finish_alignments:', 1)[0]
    snakefile = workflow_project / 'workflow/rules/worker.smk'
    snakefile.parent.mkdir(parents=True)
    environments = workflow_project / 'workflow/envs'
    environments.mkdir()
    shutil.copy2(root / 'workflow/envs/alignment.yaml', environments / 'alignment.yaml')
    snakefile.write_text(
        f"SCRIPTS = {str(root / 'workflow/scripts')!r}\n"
        "PYTHON = 'python'\nWORK = 'work/test'\nALIGNMENTS = 'results/test/orthogroups/alignments'\n"
        "LOG = 'logs/test'\nORTHOGROUP_ALIGNMENTS = 'orthogroups/alignments'\n"
        "SELECTED_ALIGNMENT_INPUTS = 'work/test/orthogroups/alignments/selected'\n"
        "ALIGNMENT_CACHE = 'work/test/orthogroups/alignments/cache'\n"
        "config = {'alignment': {'trimal_mode': 'gappyout'}}\n"
        "def alignment_batch_file(wc):\n    return 'work/test/orthogroups/alignments/batches/' + wc.batch + '.json'\n"
        + worker)
    result = subprocess.run(
        [snakemake, '--snakefile', str(snakefile), '--cores', '12',
         'results/test/orthogroups/alignments/batches/batch_000000.json'], cwd=workflow_project,
        text=True, capture_output=True, timeout=60)
    assert result.returncode == 0, result.stdout + result.stderr
    report = json.loads((workflow_project / 'logs/test/orthogroups/alignments/OG1.json').read_text())
    assert report['method'] == 'famsa'
    assert report['trimming']['mode'] == 'gappyout'
    assert report['sequences'] == 2 and 0 < report['sites'] <= report['trimming']['raw_sites']
