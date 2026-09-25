"""Translation reuse must preserve input identity and reject corrupt artifacts."""
import gzip
import json
import shutil
from pathlib import Path

import pytest
import yaml

import protein_cache
from common import write_tsv
from phase_config import write_profile
from sample_table import sample_table
from translate_cds import translate


@pytest.fixture
def seqkit():
    command = shutil.which('seqkit')
    if not command:
        pytest.skip('seqkit required')
    return command


@pytest.fixture
def cds(tmp_path):
    path = tmp_path / 'cds.fa.gz'
    with gzip.open(path, 'wt') as out:
        out.write('>Species_g1\nATGAAATAA\n')
    return path


def forbidden(*args, **kwargs):
    raise AssertionError('cached data must not invoke translation')


def test_reuse_without_seqkit_and_corruption_detection(tmp_path, cds, seqkit, monkeypatch):
    cache = tmp_path / 'cache'
    first, second = tmp_path / 'first.fa', tmp_path / 'second.fa'
    assert protein_cache.cached_translate(cds, first, first.with_suffix('.json'), cache, seqkit) is False
    monkeypatch.setattr(protein_cache, 'translate', forbidden)
    assert protein_cache.cached_translate(cds, second, second.with_suffix('.json'), cache, 'not-installed') is True
    assert first.read_bytes() == second.read_bytes()
    receipt = json.loads(second.with_suffix('.json').read_text())
    assert receipt['reused'] is True
    cached = next(cache.glob('*/table_1/protein.fa'))
    cached.write_text('>Species_g1\nWRONG\n')
    with pytest.raises(ValueError, match='registered file changed'):
        protein_cache.cached_translate(cds, second, second.with_suffix('.json'), cache)


def test_different_cds_or_genetic_code_cannot_reuse(tmp_path, cds, seqkit, monkeypatch):
    out, provenance, cache = tmp_path/'out.fa', tmp_path/'out.json', tmp_path/'cache'
    protein_cache.cached_translate(cds, out, provenance, cache, seqkit)
    monkeypatch.setattr(protein_cache, 'translate', forbidden)
    with pytest.raises(AssertionError, match='must not invoke'):
        protein_cache.cached_translate(cds, out, provenance, cache, table=2)
    with gzip.open(cds, 'wt') as handle:
        handle.write('>Species_g1\nATGCCCTAA\n')
    with pytest.raises(AssertionError, match='must not invoke'):
        protein_cache.cached_translate(cds, out, provenance, cache)
    assert len(list(cache.glob('*/table_*/receipt.json'))) == 1


def test_imported_translation_is_independent_and_checks_provenance(tmp_path, cds, seqkit, monkeypatch):
    old, provenance, cache = tmp_path/'old.fa', tmp_path/'old.json', tmp_path/'cache'
    translate(str(cds), old, provenance, seqkit)
    info = json.loads(provenance.read_text())
    info['cds']['sha256'] = '0' * 64
    wrong = tmp_path/'wrong.json'; wrong.write_text(json.dumps(info))
    with pytest.raises(ValueError, match='provenance differs'):
        protein_cache.register_translation(cds, old, wrong, cache)
    protein_cache.register_translation(cds, old, provenance, cache)
    old.unlink(); provenance.unlink()
    monkeypatch.setattr(protein_cache, 'translate', forbidden)
    assert protein_cache.cached_translate(cds, tmp_path/'new.fa', tmp_path/'new.json', cache)


def test_sample_index_refreshes_after_checkpoint_replacement(tmp_path):
    path = tmp_path/'samples.tsv'
    fields = ['odb_species', 'run']
    write_tsv(path, fields, [{'odb_species':'A', 'run':'A1'}, {'odb_species':'A', 'run':'A2'}])
    rows, species, runs = sample_table(path)
    assert len(species['A']) == 2 and runs['A2'][0]['odb_species'] == 'A'
    replacement = tmp_path/'replacement.tsv'
    write_tsv(replacement, fields, [{'odb_species':'B', 'run':'B1'}])
    replacement.replace(path)
    rows, species, runs = sample_table(path)
    assert list(species) == ['B'] and list(runs) == ['B1']


def test_slurm_profile_uses_short_squeue_polling(tmp_path):
    cfg = yaml.safe_load((Path(__file__).resolve().parents[1]/'config/build.yaml').read_text())
    profile = write_profile(tmp_path/'profile', cfg['slurm'])
    data = yaml.safe_load((profile/'config.yaml').read_text())
    assert data['jobs'] == 64
    assert 'group-components' not in data  # Snakemake 9.8 groups fail with containerized Conda.
    assert data['slurm-status-command'] == 'squeue'
    assert data['slurm-init-seconds-before-status-checks'] == 10
    assert data['seconds-between-status-checks'] == 10


def test_concurrent_builds_publish_translation_once(tmp_path, cds, seqkit, monkeypatch):
    from concurrent.futures import ThreadPoolExecutor
    original = protein_cache.translate
    calls = []
    def observed(*args, **kwargs):
        calls.append(1)
        return original(*args, **kwargs)
    monkeypatch.setattr(protein_cache, 'translate', observed)
    def run(index):
        return protein_cache.cached_translate(cds, tmp_path/f'{index}.fa', tmp_path/f'{index}.json',
                                              tmp_path/'cache', seqkit)
    with ThreadPoolExecutor(max_workers=2) as workers:
        reused = list(workers.map(run, [1, 2]))
    assert sorted(reused) == [False, True]
    assert len(calls) == 1
    assert (tmp_path/'1.fa').read_bytes() == (tmp_path/'2.fa').read_bytes()
