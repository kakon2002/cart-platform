"""Shared cache and retrieval layer.

Three properties the rest of the pipeline relies on:

* **The manifest is the commit marker.** The payload is written and moved into
  place first; the manifest is written last. A run interrupted at any point
  therefore reads as missing rather than as complete, which is the safe
  direction — a truncated file that looks finished is undetectable later.
* **The fingerprint covers the request, not just the file.** Query terms, field
  lists and release pins all feed it, so changing any of them invalidates the
  entry instead of silently serving data fetched under different terms.
* **Declared counts are checked against observed counts.** A dropped page in a
  paginated download is byte-for-byte indistinguishable from a complete one once
  it is on disk, so the check has to happen at fetch time.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import os
import time
import urllib.error
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Iterable, Iterator

CACHE_ROOT = Path(__file__).resolve().parents[2] / "data"

CHUNK = 1 << 20
MANIFEST_VERSION = 1
USER_AGENT = "cart-platform/1.0 (research pipeline)"


class CacheError(RuntimeError):
    """Base for cache and retrieval failures."""


class IntegrityError(CacheError):
    """A download completed but does not match what the source declared."""


@dataclass(frozen=True)
class CacheEntry:
    """One cached artifact and the request terms that define it."""

    key: str
    filename: str
    fingerprint: dict = field(default_factory=dict)

    def fingerprint_hash(self) -> str:
        canonical = json.dumps(self.fingerprint, sort_keys=True, separators=(",", ":"))
        return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def _atomic_replace(tmp: Path, dest: Path) -> None:
    """Move a fully written temporary file into place in one step."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    # Flush through to the device before the rename, so the payload is durable
    # by the time the manifest that blesses it can be written. The handle has to
    # be writable for this to be permitted.
    fd = os.open(tmp, os.O_RDWR | getattr(os, "O_BINARY", 0))
    try:
        os.fsync(fd)
    finally:
        os.close(fd)
    os.replace(tmp, dest)


def _write_json_atomic(dest: Path, payload: dict) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".partial")
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, sort_keys=True)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, dest)


class SourceCache:
    """Namespaced cache directory with manifest-gated reads."""

    def __init__(self, namespace: str, root: Path | None = None) -> None:
        self.namespace = namespace
        self.root = (root or CACHE_ROOT) / namespace
        self.root.mkdir(parents=True, exist_ok=True)

    def path(self, entry: CacheEntry) -> Path:
        return self.root / entry.filename

    def manifest_path(self, entry: CacheEntry) -> Path:
        return self.root / f"{entry.key}.manifest.json"

    def read_manifest(self, entry: CacheEntry) -> dict | None:
        mp = self.manifest_path(entry)
        if not mp.exists():
            return None
        try:
            with open(mp, "r", encoding="utf-8") as fh:
                return json.load(fh)
        except (OSError, json.JSONDecodeError):
            return None

    def is_valid(self, entry: CacheEntry) -> bool:
        """True only if the manifest exists, matches, and the payload is intact."""
        manifest = self.read_manifest(entry)
        if manifest is None:
            return False
        if manifest.get("manifest_version") != MANIFEST_VERSION:
            return False
        if manifest.get("fingerprint_hash") != entry.fingerprint_hash():
            return False
        payload = self.path(entry)
        if not payload.exists():
            return False
        return payload.stat().st_size == manifest.get("bytes")

    def commit(
        self,
        entry: CacheEntry,
        tmp_path: Path,
        *,
        declared_rows: int | None = None,
        observed_rows: int | None = None,
        digest: str | None = None,
        extra: dict | None = None,
    ) -> Path:
        """Move the payload into place, then write the manifest that blesses it."""
        if (
            declared_rows is not None
            and observed_rows is not None
            and declared_rows != observed_rows
        ):
            tmp_path.unlink(missing_ok=True)
            raise IntegrityError(
                f"{entry.key}: source declared {declared_rows} rows, "
                f"received {observed_rows}"
            )

        dest = self.path(entry)
        _atomic_replace(tmp_path, dest)

        manifest = {
            "manifest_version": MANIFEST_VERSION,
            "key": entry.key,
            "filename": entry.filename,
            "fingerprint": entry.fingerprint,
            "fingerprint_hash": entry.fingerprint_hash(),
            "bytes": dest.stat().st_size,
            "sha256": digest,
            "declared_rows": declared_rows,
            "observed_rows": observed_rows,
            "retrieved_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        }
        if extra:
            manifest.update(extra)
        _write_json_atomic(self.manifest_path(entry), manifest)
        return dest

    def tmp_path(self, entry: CacheEntry) -> Path:
        return self.root / f"{entry.filename}.partial"

    def ensure(
        self,
        entry: CacheEntry,
        fetcher: Callable[[Path], dict],
        *,
        force: bool = False,
    ) -> Path:
        """Return the cached payload, fetching it first if the entry is not valid.

        ``fetcher`` writes into the temporary path it is handed and returns
        manifest metadata; it is never asked to move anything into place.
        """
        if not force and self.is_valid(entry):
            return self.path(entry)

        tmp = self.tmp_path(entry)
        tmp.unlink(missing_ok=True)
        meta = fetcher(tmp) or {}
        return self.commit(entry, tmp, **meta)


# --------------------------------------------------------------------------
# retrieval
# --------------------------------------------------------------------------


def _open(
    url: str,
    *,
    method: str = "GET",
    data: bytes | None = None,
    headers: dict | None = None,
    timeout: int = 300,
    retries: int = 4,
):
    hdrs = {"User-Agent": USER_AGENT}
    if headers:
        hdrs.update(headers)
    last: Exception | None = None
    for attempt in range(retries):
        req = urllib.request.Request(url, data=data, headers=hdrs, method=method)
        try:
            return urllib.request.urlopen(req, timeout=timeout)
        except urllib.error.HTTPError as exc:
            if exc.code in (429, 500, 502, 503, 504) and attempt < retries - 1:
                last = exc
                time.sleep(2**attempt * 3)
                continue
            raise
        except (urllib.error.URLError, TimeoutError, ConnectionError) as exc:
            last = exc
            if attempt < retries - 1:
                time.sleep(2**attempt * 3)
                continue
            raise
    raise CacheError(f"unreachable: {url}") from last


def stream_to_file(
    url: str,
    dest: Path,
    *,
    method: str = "GET",
    data: bytes | None = None,
    headers: dict | None = None,
    timeout: int = 300,
    progress_label: str | None = None,
) -> dict:
    """Stream a response to disk without holding it in memory.

    Returns manifest metadata: byte count and digest.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    digest = hashlib.sha256()
    total = 0
    started = time.time()
    with _open(
        url, method=method, data=data, headers=headers, timeout=timeout
    ) as resp:
        declared_bytes = resp.headers.get("Content-Length")
        with open(dest, "wb") as fh:
            while True:
                block = resp.read(CHUNK)
                if not block:
                    break
                fh.write(block)
                digest.update(block)
                total += len(block)
                if progress_label and total % (64 * CHUNK) < CHUNK:
                    mb = total / 1_048_576
                    print(f"    {progress_label}: {mb:,.0f} MB", flush=True)
            fh.flush()
            os.fsync(fh.fileno())

    if declared_bytes is not None and int(declared_bytes) != total:
        dest.unlink(missing_ok=True)
        raise IntegrityError(
            f"{url}: declared {declared_bytes} bytes, received {total}"
        )

    if progress_label:
        secs = time.time() - started
        print(
            f"    {progress_label}: {total / 1_048_576:,.0f} MB in {secs:,.0f}s",
            flush=True,
        )
    return {"digest": digest.hexdigest(), "extra": {"source_url": url}}


def stream_paginated_to_file(
    url: str,
    dest: Path,
    *,
    headers: dict | None = None,
    timeout: int = 300,
    declared_header: str = "x-total-results",
    has_header_row: bool = True,
    progress_label: str | None = None,
) -> dict:
    """Follow ``next`` links, appending rows to one file, counting as it goes.

    The declared total is read from the first response and compared against the
    rows actually written. A dropped page is otherwise invisible.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    digest = hashlib.sha256()
    rows = 0
    declared: int | None = None
    page = 0
    next_url: str | None = url

    with open(dest, "wb") as fh:
        while next_url:
            with _open(next_url, headers=headers, timeout=timeout) as resp:
                if declared is None:
                    raw = resp.headers.get(declared_header)
                    declared = int(raw) if raw is not None else None
                body = resp.read()
                link = resp.headers.get("Link", "") or ""

            text = body.decode("utf-8")
            lines = text.splitlines(keepends=True)
            if has_header_row and lines:
                if page == 0:
                    fh.write(lines[0].encode("utf-8"))
                    digest.update(lines[0].encode("utf-8"))
                lines = lines[1:]
            for line in lines:
                fh.write(line.encode("utf-8"))
                digest.update(line.encode("utf-8"))
            rows += len(lines)

            page += 1
            if progress_label and page % 10 == 0:
                print(f"    {progress_label}: {rows:,} rows", flush=True)

            next_url = None
            for part in link.split(","):
                if 'rel="next"' in part:
                    next_url = part.split(";")[0].strip().strip("<>")
                    break

        fh.flush()
        os.fsync(fh.fileno())

    if progress_label:
        print(f"    {progress_label}: {rows:,} rows over {page} pages", flush=True)

    return {
        "digest": digest.hexdigest(),
        "declared_rows": declared,
        "observed_rows": rows,
        "extra": {"source_url": url, "pages": page},
    }


def decompress_gzip(src: Path, dest: Path, *, progress_label: str | None = None) -> dict:
    """Expand a gzip payload in fixed-size blocks.

    Used for archives far larger than memory, so neither side is ever fully
    resident.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    digest = hashlib.sha256()
    total = 0
    with gzip.open(src, "rb") as fin, open(dest, "wb") as fout:
        while True:
            block = fin.read(CHUNK)
            if not block:
                break
            fout.write(block)
            digest.update(block)
            total += len(block)
            if progress_label and total % (256 * CHUNK) < CHUNK:
                print(f"    {progress_label}: {total / 1_048_576:,.0f} MB", flush=True)
        fout.flush()
        os.fsync(fout.fileno())
    if progress_label:
        print(f"    {progress_label}: {total / 1_048_576:,.0f} MB expanded", flush=True)
    return {"digest": digest.hexdigest(), "extra": {"expanded_bytes": total}}


def extract_from_zip(archive: Path, member: str, dest: Path) -> dict:
    """Pull one member out of a zip archive in blocks."""
    import zipfile

    dest.parent.mkdir(parents=True, exist_ok=True)
    digest = hashlib.sha256()
    total = 0
    with zipfile.ZipFile(archive) as zf:
        names = zf.namelist()
        if member not in names:
            matches = [n for n in names if n.endswith(member)]
            if not matches:
                raise IntegrityError(f"{member} not present in {archive.name}")
            member = matches[0]
        with zf.open(member) as fin, open(dest, "wb") as fout:
            while True:
                block = fin.read(CHUNK)
                if not block:
                    break
                fout.write(block)
                digest.update(block)
                total += len(block)
            fout.flush()
            os.fsync(fout.fileno())
    return {"digest": digest.hexdigest(), "extra": {"member": member}}


def post_json(url: str, payload: dict, *, timeout: int = 300) -> bytes:
    """Small POST whose response fits in memory comfortably."""
    body = json.dumps(payload).encode("utf-8")
    with _open(
        url,
        method="POST",
        data=body,
        headers={"Content-Type": "application/json"},
        timeout=timeout,
    ) as resp:
        return resp.read()


def read_lines(path: Path, encoding: str = "utf-8") -> Iterator[str]:
    """Iterate a cached text payload one line at a time."""
    with open(path, "r", encoding=encoding, newline="") as fh:
        for line in fh:
            yield line.rstrip("\r\n")


# --------------------------------------------------------------------------
# source contract
# --------------------------------------------------------------------------


class DataSource:
    """Base class. A source is usable only when every entry it declares is valid."""

    name: str = "unnamed"
    namespace: str = "unnamed"

    def __init__(self, root: Path | None = None) -> None:
        self.cache = SourceCache(self.namespace, root=root)

    def cache_entries(self) -> Iterable[CacheEntry]:
        raise NotImplementedError

    def is_cached(self) -> bool:
        """A partial cache is not a usable source."""
        entries = list(self.cache_entries())
        return bool(entries) and all(self.cache.is_valid(e) for e in entries)

    def missing_entries(self) -> list[CacheEntry]:
        return [e for e in self.cache_entries() if not self.cache.is_valid(e)]

    def fetch(self) -> None:
        raise NotImplementedError
