"""Optional checksum memoization, guarded by the complete filesystem identity."""
from contextlib import contextmanager
from contextvars import ContextVar
import json
import os
from pathlib import Path
import sqlite3


_ACTIVE = ContextVar("verification_cache", default=None)


class VerificationCache:
    def __init__(self, path):
        self.files, self.pending = {}, {}
        self.connection = None
        try:
            path = Path(path)
            path.parent.mkdir(parents=True, exist_ok=True)
            self.connection = sqlite3.connect(path, timeout=1)
            self.connection.execute("CREATE TABLE IF NOT EXISTS verified "
                                    "(path TEXT, sha256 TEXT, stat TEXT, PRIMARY KEY(path, sha256))")
            self.connection.commit()
            for filename, checksum, identity in self.connection.execute("SELECT path, sha256, stat FROM verified"):
                try:
                    self.files[filename, checksum] = tuple(json.loads(identity))
                except (ValueError, TypeError):
                    continue
        except (OSError, sqlite3.Error):
            self._disconnect()

    def _disconnect(self):
        if self.connection is not None:
            self.connection.close()
            self.connection = None

    def matches(self, path, checksum, identity):
        return self.files.get((os.path.abspath(path), checksum)) == tuple(identity)

    def remember(self, path, checksum, identity):
        key = (os.path.abspath(path), checksum)
        self.files[key] = tuple(identity)
        self.pending[key] = tuple(identity)
        if len(self.pending) >= 1000:
            self.flush()

    def flush(self):
        if self.connection is not None and self.pending:
            try:
                # Transactions only cover this short batch; never hold a write lock while hashing.
                with self.connection:
                    self.connection.executemany("INSERT OR REPLACE INTO verified VALUES (?, ?, ?)",
                        [(path, checksum, json.dumps(identity))
                         for (path, checksum), identity in self.pending.items()])
            except sqlite3.Error:
                self._disconnect()
        self.pending.clear()

    def close(self):
        self.flush()
        self._disconnect()


def active_cache():
    return _ACTIVE.get()


@contextmanager
def verification_cache(path):
    """Share checks within a command and persist successful hashes across commands."""
    if _ACTIVE.get() is not None:
        yield _ACTIVE.get()
        return
    cache = VerificationCache(path)
    token = _ACTIVE.set(cache)
    try:
        yield cache
    finally:
        _ACTIVE.reset(token)
        cache.close()
