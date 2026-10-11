"""Publication carries verified identities through hard links without rereading payloads."""
import gzip
import json
from pathlib import Path

import pytest

import common
import dataset
import dataset_assets as assets
from product_storage import product_record
from verification_cache import verification_cache


@pytest.mark.parametrize('compressed_busco', [False, True])
def test_materialize_reuses_verified_bytes_and_records_final_paths(tmp_path, monkeypatch, compressed_busco):
    source = tmp_path / 'source'
    source.mkdir()
    cds, abundance = source / 'cds.fa.gz', source / 'abundance.tsv.gz'
    for path in (cds, abundance):
        with gzip.open(path, 'wb') as handle:
            handle.write(b'verified payload\n')
    busco = source / ('busco.tsv.gz' if compressed_busco else 'busco.tsv')
    if compressed_busco:
        with gzip.open(busco, 'wb') as handle:
            handle.write(b'full BUSCO table\n')
    else:
        busco.write_bytes(b'full BUSCO table\n')
    entries = [product_record(path) for path in (cds, abundance, busco)]
    products = {'reference': {'cds': entries[0]}, 'quant': {'abundance': entries[1]},
                'busco': {'full': entries[2], 'counts': dict.fromkeys(assets.COUNTS, 1)}}
    item = {'species': 'Plant_A', 'row': {'species': 'Plant_A', 'run': 'A1'}}
    target = tmp_path / 'input'
    manifest = {'items': [item], 'analysis': {'input_root': str(target)}, 'fields': ['species', 'run']}
    monkeypatch.setattr(dataset, 'load', lambda *args, **kwargs: manifest)
    monkeypatch.setattr(dataset, 'item_products', lambda *args: products)
    identities = {tuple(entry['stat'][:2]) for entry in entries}
    hash_file = common.sha256
    def forbid_payload_reads(path):
        stat = Path(path).stat()
        if (stat.st_dev, stat.st_ino) in identities:
            pytest.fail(f'verified hard-linked payload was rehashed: {path}')
        return hash_file(path)
    monkeypatch.setattr(common, 'sha256', forbid_payload_reads)
    monkeypatch.setattr(assets, 'sha256', forbid_payload_reads)
    with verification_cache(tmp_path / 'verification.sqlite') as cache:
        for entry in entries:
            cache.remember(entry['path'], entry['sha256'], entry['stat'])
        assert dataset.materialize(tmp_path) == target
        for entry in entries:
            assets.verify(entry)
        assert dataset.materialize(tmp_path) == target
    assert cds.samefile(target / 'cds/Plant_A_longestCDS.fa.gz')
    assert abundance.samefile(target / 'quant/Plant_A/A1/A1_abundance.tsv.gz')
    assert (target / 'busco/full/Plant_A.busco.full.tsv').read_bytes() == b'full BUSCO table\n'
    receipt = json.loads((tmp_path / 'input_receipt.json').read_text())
    assert len(receipt['files']) == 5
    assert all(Path(entry['path']).is_relative_to(target) for entry in receipt['files'])
    assert all(entry['stat'] == assets.stat_identity(entry['path']) for entry in receipt['files'])
    assert entries[1]['content_sha256'] == next(
        entry['content_sha256'] for entry in receipt['files'] if '/quant/' in entry['path'])
