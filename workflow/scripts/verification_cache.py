"""Optional checksum and input indexes, guarded by the complete file identity."""
from contextlib import contextmanager
from contextvars import ContextVar
import hashlib
import json
import os
from pathlib import Path
import sqlite3
import zlib


_ACTIVE = ContextVar("verification_cache", default=None)


class VerificationCache:
    def __init__(self, path):
        self.files, self.pending, self.identities, self.checksums = {}, {}, set(), {}
        self.connection = None
        try:
            path = Path(path)
            path.parent.mkdir(parents=True, exist_ok=True)
            self.connection = sqlite3.connect(path, timeout=1)
            self.connection.execute("CREATE TABLE IF NOT EXISTS verified "
                                    "(path TEXT, sha256 TEXT, stat TEXT, PRIMARY KEY(path, sha256))")
            self.connection.execute("CREATE TABLE IF NOT EXISTS documents "
                                    "(path TEXT, kind TEXT, source_sha256 TEXT, stat TEXT, "
                                    "payload_sha256 TEXT, payload BLOB, PRIMARY KEY(path, kind))")
            self.connection.commit()
            for filename, checksum, identity in self.connection.execute("SELECT path, sha256, stat FROM verified"):
                try:
                    identity = tuple(json.loads(identity))
                    self.files[filename, checksum] = identity
                    self.identities.add((checksum, identity))
                    self.checksums[identity] = checksum
                except (ValueError, TypeError):
                    continue
        except (OSError, sqlite3.Error):
            self._disconnect()

    def _disconnect(self):
        if self.connection is not None:
            self.connection.close()
            self.connection = None

    def matches(self, path, checksum, identity):
        # Hard-link aliases are the same bytes only when the entire identity,
        # including ctime, agrees. Never reuse a match across a metadata change.
        return (checksum, tuple(identity)) in self.identities

    def remember(self, path, checksum, identity):
        key = (os.path.abspath(path), checksum)
        if self.files.get(key) == tuple(identity):
            return
        self.files[key] = tuple(identity)
        self.identities.add((checksum, tuple(identity)))
        self.checksums[tuple(identity)] = checksum
        self.pending[key] = tuple(identity)
        if len(self.pending) >= 1000:
            self.flush()

    def checksum(self, identity):
        return self.checksums.get(tuple(identity))

    def document(self, path, kind, identity, checksum=None):
        """Read a versioned index only while its source and payload still agree."""
        if self.connection is None:
            return None
        try:
            row = self.connection.execute(
                "SELECT source_sha256, stat, payload_sha256, payload FROM documents "
                "WHERE path = ? AND kind = ?", (os.path.abspath(path), kind)).fetchone()
            if row is None:
                return None
            source, saved, expected, payload = row
            if tuple(json.loads(saved)) != tuple(identity) or (checksum is not None and source != checksum):
                return None
            if hashlib.sha256(payload).hexdigest() != expected:
                return None
            data = json.loads(zlib.decompress(payload))
            return data if isinstance(data, dict) else None
        except sqlite3.Error:
            self._disconnect()
        except (ValueError, TypeError, zlib.error):
            pass  # A damaged optional index is rebuilt from its source.
        return None

    def remember_document(self, path, kind, identity, checksum, data):
        if self.connection is None:
            return
        # Serialization happens before the short transaction, outside any write lock.
        payload = zlib.compress(json.dumps(data, separators=(",", ":")).encode(), level=1)
        expected = hashlib.sha256(payload).hexdigest()
        try:
            with self.connection:
                self.connection.execute("INSERT OR REPLACE INTO documents VALUES (?, ?, ?, ?, ?, ?)",
                    (os.path.abspath(path), kind, checksum, json.dumps(identity), expected, payload))
        except sqlite3.Error:
            self._disconnect()

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


def compact_rows(rows):
    """Share column definitions; omit only empty strings, preserving missing keys."""
    columns, schemas, packed = [], {}, []
    for row in rows:
        keys = tuple(row)
        index = schemas.get(keys)
        if index is None:
            index = len(columns)
            schemas[keys] = index
            columns.append(list(keys))
        packed.append({'columns': index, 'values': {k: v for k, v in row.items() if v != ''}})
    return columns, packed


def expand_rows(columns, rows):
    """Reconstruct the original row keys and empty values for existing consumers."""
    for row in rows:
        expanded = dict.fromkeys(columns[row['columns']], '')
        expanded.update(row['values'])
        yield expanded


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
