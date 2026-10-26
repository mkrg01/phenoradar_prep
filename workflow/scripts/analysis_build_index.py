"""Small, verified product indexes for downstream workers (no metadata rows)."""
import json
from collections.abc import Mapping
from pathlib import Path

from common import write_json
from dataset_assets import digest, record, verify


def create_index(completion, destination, completed=None, species=None):
    from build_products import load_complete
    from portable_build import completion_path
    source = Path(completion).resolve() if completed is not None else completion_path(completion)
    completed = completed if completed is not None else load_complete(source, verify_files=False)
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    products = {}
    selected = set(completed['products']) if species is None else set(species)
    if selected - completed['products'].keys():
        raise ValueError('analysis index selection is absent from the database')
    for number, species in enumerate(sorted(selected)):
        product = completed['products'][species]
        # Importers need the run identity and file receipts, not the metadata row.
        small = {k: product[k] for k in ('odb_species', 'cds', 'abundance', 'protein',
                                       'translation', 'expression', 'expression_qc')}
        small['row'] = {'run': product['row']['run']}
        path = destination / 'products' / f'{number:06d}.json'
        write_json(path, small)
        entry = record(path)
        entry['path'] = str(path.relative_to(destination))
        products[species] = {'odb_species': product['odb_species'],
                             'run': product['row']['run'], 'file': entry}
    index = {'schema_version': 1, 'source': record(source), 'mapping': completed['mapping'],
             'tpm': completed['tpm'], 'products': products}
    index['sha256'] = digest(index)
    write_json(destination / 'index.json', index)
    return destination / 'index.json'


class Products(Mapping):
    def __init__(self, root, entries):
        self.root, self.entries = root, entries

    def __iter__(self):
        return iter(self.entries)

    def __len__(self):
        return len(self.entries)

    def __getitem__(self, species):
        entry = self.entries[species]['file']
        relative = Path(entry['path'])
        if relative.is_absolute() or '..' in relative.parts:
            raise ValueError('invalid analysis product index path')
        path = verify(dict(entry, path=str(self.root / relative)))
        product = json.loads(path.read_text())
        identity = self.entries[species]
        if product['odb_species'] != identity['odb_species'] or product['row']['run'] != identity['run']:
            raise ValueError('analysis product index identity changed')
        return product


def load_index(path, completion):
    path = Path(path)
    data = json.loads(path.read_text())
    if data.get('schema_version') != 1 or digest({k: v for k, v in data.items() if k != 'sha256'}) != data.get('sha256'):
        raise ValueError('analysis build index changed')
    # Unchanged stat receipts avoid hashing/parsing the large source on workers.
    if Path(data['source']['path']).resolve() != Path(completion).resolve():
        raise ValueError('analysis index belongs to another database')
    verify(data['source'])
    return dict(data, products=Products(path.parent, data['products']))


class BuildInputs:
    """Resolve imports lazily; unrelated workers never open the database/index."""
    def __init__(self, completion, index):
        self.completion, self.index, self._data = completion, Path(index), None

    @property
    def data(self):
        if self._data is None:
            from build_products import load_complete
            self._data = (load_index(self.index, self.completion) if self.index.exists()
                          else load_complete(self.completion, verify_files=False))
        return self._data

    def product(self, value, field):
        products = self.data['products']
        if isinstance(products, Products):
            matches = [s for s, p in products.entries.items() if p[field] == value]
        else:
            matches = [s for s, p in products.items()
                       if (p['row']['run'] if field == 'run' else p[field]) == value]
        if len(matches) != 1:
            raise ValueError(f'expected one database product for {field}={value}')
        return products[matches[0]]
