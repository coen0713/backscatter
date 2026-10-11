"""Resumable downloads with retries.

Large public files (DEM tiles, the geoid grid, WorldCover) are fetched once into
data/cache. Connections that are reset mid-transfer are resumed with an HTTP
Range request, so a flaky network costs time, not correctness.
"""

from __future__ import annotations

import time
import urllib.error
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CACHE = ROOT / "data" / "cache"


def _content_length(url: str, retries: int) -> int | None:
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, method="HEAD")
            with urllib.request.urlopen(req, timeout=60) as resp:
                n = resp.headers.get("Content-Length")
                return int(n) if n is not None else None
        except (urllib.error.URLError, ConnectionError, TimeoutError):
            time.sleep(min(2.0 ** attempt, 30.0))
    return None


def download(url: str, dest: Path | None = None, retries: int = 12) -> Path:
    """Download `url` to `dest` (default: data/cache/<basename>) unless it is
    already complete. Returns the local path."""
    dest = dest or CACHE / url.rsplit("/", 1)[1]
    dest.parent.mkdir(parents=True, exist_ok=True)
    total = _content_length(url, retries)
    if dest.exists() and (total is None or dest.stat().st_size == total):
        return dest
    part = dest.with_suffix(dest.suffix + ".part")
    for attempt in range(retries):
        have = part.stat().st_size if part.exists() else 0
        if total is not None and have >= total:
            break
        headers = {"Range": f"bytes={have}-"} if have else {}
        try:
            req = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(req, timeout=60) as resp, open(part, "ab" if have else "wb") as f:
                if have and resp.status != 206:  # server ignored the range: start over
                    f.truncate(0)
                while chunk := resp.read(1 << 20):
                    f.write(chunk)
            if total is None or part.stat().st_size >= total:
                break
        except (urllib.error.URLError, ConnectionError, TimeoutError) as err:
            print(f"  retry {attempt + 1}/{retries} for {dest.name}: {err}", flush=True)
            time.sleep(min(2.0 ** attempt, 30.0))
    else:
        raise RuntimeError(f"download failed after {retries} attempts: {url}")
    if total is not None and part.stat().st_size != total:
        raise RuntimeError(f"incomplete download of {url}: {part.stat().st_size} of {total} bytes")
    part.replace(dest)
    return dest
