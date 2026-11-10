"""Successful expression checks are reusable only for the same full contract."""
import hashlib
import json
import os
import shutil
import sqlite3
import zlib

import pytest

import analysis_validation as validation
from common import write_tsv
from verification_cache import verification_cache


RUNS = {'R1': 'Plant_A', 'R2': 'Plant_B'}


@pytest.fixture
def expression(tmp_path):
    path = tmp_path / 'tpm.tsv'
    write_tsv(path, ['species', 'orthogroup', 'tpm'], [
        dict(species=sample, orthogroup='OG1', tpm='1.0') for sample in RUNS.values()])
    return path


def check(path, **kwargs):
    return validation.validate_expression(path, RUNS, 'orthogroup', 'tpm', **kwargs)


def test_successful_check_persists_without_reopening_rows(tmp_path, expression, monkeypatch):
    cache = tmp_path / 'cache.sqlite'
    with verification_cache(cache):
        assert check(expression) == {'OG1'}
    monkeypatch.setattr(validation, '_validate_expression', lambda *a, **k: pytest.fail('reread rows'))
    with verification_cache(cache):
        result = check(expression)
        result.add('caller mutation')
        assert check(expression) == {'OG1'}


@pytest.mark.parametrize('storage', ['copy', 'hardlink'])
def test_verified_content_reuses_check_at_new_path(tmp_path, expression, monkeypatch, storage):
    cache = tmp_path / 'cache.sqlite'
    with verification_cache(cache):
        check(expression)
    alias = tmp_path / 'other.tsv'
    (shutil.copy2 if storage == 'copy' else os.link)(expression, alias)
    monkeypatch.setattr(validation, '_validate_expression', lambda *a, **k: pytest.fail('reread rows'))
    with verification_cache(cache):
        assert check(alias) == {'OG1'}


def test_restored_size_and_mtime_do_not_hide_bad_values(tmp_path, expression):
    cache = tmp_path / 'cache.sqlite'
    with verification_cache(cache):
        check(expression)
    before = expression.stat()
    expression.write_text(expression.read_text().replace('1.0', 'NaN'))
    os.utime(expression, ns=(before.st_atime_ns, before.st_mtime_ns))
    with verification_cache(cache), pytest.raises(ValueError, match='finite and nonnegative'):
        check(expression)


@pytest.mark.parametrize('change', ['runs', 'feature', 'value', 'validator', 'collect_features'])
def test_changed_contract_requires_row_validation(tmp_path, expression, monkeypatch, change):
    cache = tmp_path / 'cache.sqlite'
    with verification_cache(cache):
        assert check(expression, collect_features=False) == set()
    original = validation._validate_expression
    calls = []
    def counted(*args, **kwargs):
        calls.append(True)
        return original(*args, **kwargs)
    monkeypatch.setattr(validation, '_validate_expression', counted)
    runs, feature, value, collect = dict(RUNS), 'orthogroup', 'tpm', False
    if change == 'runs':
        runs['R3'] = runs.pop('R1')
    elif change == 'feature':
        feature = 'other_feature'
    elif change == 'value':
        value = 'other_value'
    elif change == 'validator':
        monkeypatch.setattr(validation, '_EXPRESSION_VALIDATOR_SHA', 'new validator')
    else:
        collect = True
    with verification_cache(cache):
        if change in {'feature', 'value'}:
            with pytest.raises(ValueError, match='missing/duplicate columns'):
                validation.validate_expression(expression, runs, feature, value, collect_features=collect)
        else:
            result = validation.validate_expression(expression, runs, feature, value, collect_features=collect)
            assert result == ({'OG1'} if collect else set())
    assert len(calls) == 1


def test_compression_mode_is_part_of_contract(tmp_path, expression):
    cache = tmp_path / 'cache.sqlite'
    with verification_cache(cache):
        check(expression)
    alias = tmp_path / 'other.tsv.gz'
    shutil.copy2(expression, alias)
    with verification_cache(cache), pytest.raises(OSError):
        check(alias)


def test_failed_validation_is_never_cached(tmp_path, expression, monkeypatch):
    expression.write_text(expression.read_text().replace('1.0', 'NaN'))
    original, calls = validation._validate_expression, []
    def counted(*args, **kwargs):
        calls.append(True)
        return original(*args, **kwargs)
    monkeypatch.setattr(validation, '_validate_expression', counted)
    for _ in range(2):
        with verification_cache(tmp_path / 'cache.sqlite'), pytest.raises(ValueError):
            check(expression)
    assert len(calls) == 2


def test_changes_during_validation_do_not_publish_success(tmp_path, expression, monkeypatch):
    original = validation._validate_expression
    def changing(*args, **kwargs):
        result = original(*args, **kwargs)
        expression.write_text(expression.read_text().replace('1.0', 'NaN'))
        return result
    monkeypatch.setattr(validation, '_validate_expression', changing)
    with verification_cache(tmp_path / 'cache.sqlite'), pytest.raises(ValueError, match='changed during'):
        check(expression)
    monkeypatch.setattr(validation, '_validate_expression', original)
    with verification_cache(tmp_path / 'cache.sqlite'), pytest.raises(ValueError, match='finite and nonnegative'):
        check(expression)


def test_cache_hit_still_checks_source_after_lookup(tmp_path, expression, monkeypatch):
    with verification_cache(tmp_path / 'cache.sqlite') as cache:
        check(expression)
        original = cache.validation
        def changing(*args):
            result = original(*args)
            expression.write_text(expression.read_text().replace('1.0', 'NaN'))
            return result
        monkeypatch.setattr(cache, 'validation', changing)
        with pytest.raises(ValueError, match='changed during'):
            check(expression)


@pytest.mark.parametrize('damage', ['checksum', 'compressed_data', 'schema'])
def test_damaged_optional_proof_is_rebuilt(tmp_path, expression, monkeypatch, damage):
    cache = tmp_path / 'cache.sqlite'
    with verification_cache(cache):
        check(expression)
    payload = zlib.compress(json.dumps({'features': None}).encode()) if damage == 'schema' else b'broken'
    checksum = '0' * 64 if damage == 'checksum' else hashlib.sha256(payload).hexdigest()
    with sqlite3.connect(cache) as connection:
        connection.execute('UPDATE validations SET payload = ?, payload_sha256 = ?', (payload, checksum))
    calls, original = [], validation._validate_expression
    def counted(*args, **kwargs):
        calls.append(True)
        return original(*args, **kwargs)
    monkeypatch.setattr(validation, '_validate_expression', counted)
    with verification_cache(cache):
        assert check(expression) == {'OG1'}
        assert check(expression) == {'OG1'}
    assert len(calls) == 1


def test_unavailable_cache_preserves_validation(tmp_path, expression):
    cache = tmp_path / 'broken.sqlite'
    cache.write_bytes(b'not a database')
    with verification_cache(cache):
        assert check(expression) == {'OG1'}
        expression.write_text(expression.read_text().replace('1.0', 'NaN'))
        with pytest.raises(ValueError, match='finite and nonnegative'):
            check(expression)
