"""OrthoDB release selection and the shared reference pin for one build."""
import fcntl
import json
import re
from pathlib import Path

from common import write_json


def validate(odb):
    version = odb.get('version', 'v12')
    if not isinstance(version, str) or not re.fullmatch(r'v[1-9][0-9]*', version):
        raise ValueError('odb.version must be a fixed release such as v12')


def identity(odb):
    return {'version': odb.get('version', 'v12'), 'node': odb['ncbi_tax_id']}


def matches(source, odb):
    return (source.get('version'), source.get('node')) == (
        odb.get('version', 'v12'), odb['ncbi_tax_id'])


def cache_namespace(odb):
    return f'{odb.get("version", "v12")}_{odb["ncbi_tax_id"]}'


def reference_directory(odb):
    return 'resources/orthodb/' + cache_namespace(odb)


def pin_reference(path, version, node, reference_sha256=None):
    """All concurrent chunks and retries must use the same reference bytes."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(str(path) + '.lock', 'a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        expected = {'version': version, 'node': node, 'reference_sha256': reference_sha256}
        if path.exists():
            saved = json.loads(path.read_text())
            if (saved['version'], saved['node']) != (version, node):
                raise ValueError('OrthoDB reference pin version/node changed')
            previous = saved.get('reference_sha256')
            if previous and reference_sha256 and previous != reference_sha256:
                raise ValueError('OrthoDB reference differs from the pinned build snapshot; restore the original reference')
            if previous or reference_sha256 is None:
                return saved
        write_json(path, expected)
        return expected
