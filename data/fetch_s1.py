#!/usr/bin/env python3
"""Download Sentinel-1 products from the Copernicus Data Space Ecosystem.

Run this yourself, in your own terminal: it needs your CDSE account. The
username and password are read from the environment variables CDSE_USERNAME and
CDSE_PASSWORD, or asked for interactively (the password is not echoed). They
are only sent to the CDSE login service; the script never prints, logs or
writes them, and keeps the short-lived access token in memory only.

Downloads are resumable (re-run after an interruption) and verified against
the catalogue's size. Do NOT open the downloaded scenes before docs/prereg.md is
registered: they cover held-out tiles.

PowerShell:
  $env:CDSE_USERNAME = "you@example.com"
  python data/fetch_s1.py --check                  # login test only
  python data/fetch_s1.py S1A_IW_GRDH_1SDV_..._4650.SAFE S1A_IW_GRDH_1SDV_..._4683.SAFE
  python data/fetch_s1.py --selected build/validation/selected_scenes.json
"""

from __future__ import annotations

import argparse
import getpass
import json
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

TOKEN_URL = "https://identity.dataspace.copernicus.eu/auth/realms/CDSE/protocol/openid-connect/token"
CATALOGUE = "https://catalogue.dataspace.copernicus.eu/odata/v1/Products"
DOWNLOAD = "https://zipper.dataspace.copernicus.eu/odata/v1/Products({id})/$value"
ROOT = Path(__file__).resolve().parent.parent


class Session:
    """Holds the access token; refreshes it when it expires (~10 min)."""

    def __init__(self, username: str, password: str):
        self._token, self._refresh, self._expires = self._request(
            {"grant_type": "password", "username": username, "password": password,
             "client_id": "cdse-public"})

    @staticmethod
    def _request(form: dict) -> tuple[str, str, float]:
        data = urllib.parse.urlencode(form).encode()
        try:
            with urllib.request.urlopen(urllib.request.Request(TOKEN_URL, data=data), timeout=60) as r:
                body = json.load(r)
        except urllib.error.HTTPError as err:
            # Never echo the request (it contains the password).
            raise SystemExit(f"CDSE login failed (HTTP {err.code}). Check the username/password; "
                             "the account must be activated by e-mail.") from None
        return body["access_token"], body.get("refresh_token", ""), time.time() + body["expires_in"] - 30

    def token(self) -> str:
        if time.time() > self._expires:
            self._token, self._refresh, self._expires = self._request(
                {"grant_type": "refresh_token", "refresh_token": self._refresh,
                 "client_id": "cdse-public"})
        return self._token


def lookup(name: str) -> tuple[str, int]:
    """Product name -> (catalogue Id, size in bytes), via the public catalogue."""
    flt = urllib.parse.quote(f"Name eq '{name}'")
    with urllib.request.urlopen(f"{CATALOGUE}?$filter={flt}", timeout=60) as r:
        hits = json.load(r)["value"]
    if len(hits) != 1:
        raise SystemExit(f"{name}: {len(hits)} catalogue matches")
    return hits[0]["Id"], int(hits[0].get("ContentLength", 0))


def download(session: Session, name: str, out: Path, retries: int = 20) -> Path:
    pid, size = lookup(name)
    dest = out / (name.removesuffix(".SAFE") + ".zip")
    part = dest.with_suffix(".zip.part")
    if dest.exists() and (size == 0 or dest.stat().st_size == size):
        print(f"already complete: {dest.name}")
        return dest
    out.mkdir(parents=True, exist_ok=True)
    for attempt in range(retries):
        have = part.stat().st_size if part.exists() else 0
        if size and have >= size:
            break
        req = urllib.request.Request(DOWNLOAD.format(id=pid))
        req.add_header("Authorization", f"Bearer {session.token()}")
        if have:
            req.add_header("Range", f"bytes={have}-")
        try:
            with urllib.request.urlopen(req, timeout=120) as r, open(part, "ab" if have else "wb") as f:
                if have and r.status != 206:
                    f.truncate(0)
                    have = 0
                last = time.time()
                while chunk := r.read(4 << 20):
                    f.write(chunk)
                    have += len(chunk)
                    if time.time() - last > 5:
                        pct = f" ({100 * have / size:.0f}%)" if size else ""
                        print(f"  {dest.name}: {have / 1e9:.2f} GB{pct}", flush=True)
                        last = time.time()
            if not size or part.stat().st_size >= size:
                break
        except (urllib.error.URLError, ConnectionError, TimeoutError) as err:
            code = getattr(err, "code", None)
            if code in (401, 403):
                raise SystemExit(f"{name}: download refused (HTTP {code}); is the account active?")
            print(f"  retry {attempt + 1}/{retries}: {err}", flush=True)
            time.sleep(min(2.0 ** attempt, 30.0))
    else:
        raise SystemExit(f"{name}: download failed after {retries} attempts")
    if size and part.stat().st_size != size:
        raise SystemExit(f"{name}: got {part.stat().st_size} bytes, catalogue says {size}")
    part.replace(dest)
    print(f"done: {dest} ({dest.stat().st_size / 1e9:.2f} GB)")
    return dest


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("products", nargs="*", help="product names (ending in .SAFE)")
    p.add_argument("--selected", type=Path, help="JSON from validation.select_scenes --selected")
    p.add_argument("--out", type=Path, default=ROOT / "data" / "raw")
    p.add_argument("--check", action="store_true", help="only test the login")
    args = p.parse_args()

    username = os.environ.get("CDSE_USERNAME") or input("CDSE username (e-mail): ")
    password = os.environ.get("CDSE_PASSWORD") or getpass.getpass("CDSE password (not shown): ")
    session = Session(username, password)
    del password
    print("login ok")
    if args.check:
        return 0
    names = list(args.products)
    if args.selected:
        sel = json.loads(args.selected.read_text())
        names += [n for acq in sel.values() for n in acq["products"]]
    if not names:
        raise SystemExit("nothing to download: give product names or --selected")
    for name in names:
        download(session, name, args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
