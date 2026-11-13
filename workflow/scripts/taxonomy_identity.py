"""Resolve the current species and lineage in a fixed NCBI taxonomy snapshot."""
from functools import lru_cache
from pathlib import Path
import sqlite3

from dataset_assets import stat_identity


class TaxonomyIdentity:
    def __init__(self, path):
        self.path = Path(path).resolve()
        self.stat = stat_identity(self.path)
        self.db = sqlite3.connect(self.path.as_uri() + '?mode=ro', uri=True)
        self.db.execute('BEGIN')
        try:
            if self.node(1) is None:
                raise ValueError('taxonomy snapshot has no root node')
            self.db.execute('SELECT taxid_old, taxid_new FROM merged LIMIT 0')
        except Exception:
            self.db.close()
            raise

    def __enter__(self):
        return self

    def __exit__(self, exc_type, *_):
        self.db.close()
        self.node.cache_clear()
        self.resolve.cache_clear()
        if exc_type is None and stat_identity(self.path) != self.stat:
            raise ValueError('taxonomy snapshot changed during resolution')

    @lru_cache(maxsize=65536)
    def node(self, taxid):
        return self.db.execute(
            'SELECT taxid, parent, spname, rank FROM species WHERE taxid=?', (taxid,)
        ).fetchone()

    @lru_cache(maxsize=65536)
    def resolve(self, taxid):
        """Return the live taxon and named ranks, following redirects and parents."""
        taxid = int(taxid)
        seen = set()
        while True:
            if taxid in seen:
                raise ValueError('cycle in merged taxonomy IDs')
            seen.add(taxid)
            target = self.db.execute('SELECT taxid_new FROM merged WHERE taxid_old=?', (taxid,)).fetchone()
            if target is None:
                break
            taxid = target[0]
        node = self.node(taxid)
        if node is None:
            return None
        lineage = {}
        seen.clear()
        while node is not None:
            current, parent, name, rank = node
            if current in seen:
                raise ValueError(f'cycle in taxonomy lineage for {taxid}')
            seen.add(current)
            lineage.setdefault(rank, (str(current), name))
            if current == parent or current == 1:
                return str(taxid), lineage
            node = self.node(parent)
        raise ValueError(f'incomplete taxonomy lineage for {taxid}')
