#!/usr/bin/env python3
"""Download the Sentinel-1 precise orbit (POEORB) file covering a given time.

Orbit files are served without authentication from ESA's SNAP auxdata mirror:
https://step.esa.int/auxdata/orbits/Sentinel-1/POEORB/<S1A|S1B|S1C>/<YYYY>/<MM>/
Each POEORB file is valid from ~1 day before to ~1 day after its nominal day,
with state vectors every 10 s. (They are also available from the Copernicus
Data Space Ecosystem, which needs a free account.)

Example:
    python data/fetch_orbit.py --mission S1A --time 2024-06-01T05:30:00 --out data/raw/
    bsar_render --scene dem:data/raw/eiger.asc --geographic \
        --orbit data/raw/S1A_OPER_AUX_POEORB_....EOF --time 2024-06-01T05:30:00
"""

from __future__ import annotations

import argparse
import re
import sys
import time
import urllib.error
import urllib.request
import zipfile
from datetime import datetime, timedelta
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from validation.fetch import download  # noqa: E402

BASE = "https://step.esa.int/auxdata/orbits/Sentinel-1/POEORB"
PATTERN = re.compile(r"(S1[ABC]_OPER_AUX_POEORB_OPOD_\d{8}T\d{6}_V(\d{8}T\d{6})_(\d{8}T\d{6})\.EOF(?:\.zip)?)")


def listing(mission: str, month: datetime, retries: int = 8) -> list[tuple[str, datetime, datetime]]:
    url = f"{BASE}/{mission}/{month:%Y}/{month:%m}/"
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=60) as resp:
                html = resp.read().decode()
            break
        except (urllib.error.URLError, ConnectionError, TimeoutError):
            if attempt == retries - 1:
                raise
            time.sleep(min(2.0 ** attempt, 30.0))
    out = []
    for name, start, stop in sorted(set(PATTERN.findall(html))):
        out.append((url + name, datetime.strptime(start, "%Y%m%dT%H%M%S"),
                    datetime.strptime(stop, "%Y%m%dT%H%M%S")))
    return out


def fetch_orbit(mission: str, when: datetime, out: Path) -> Path:
    """Download and unzip the POEORB file covering `when` (+-10 min) into `out`;
    returns the .EOF path (reused if already present)."""
    candidates = []
    months = {(when + timedelta(days=d)).replace(day=1, hour=0, minute=0, second=0, microsecond=0)
              for d in (-2, 0, 2)}
    for month in sorted(months):
        try:
            candidates += listing(mission, month)
        except OSError as err:
            print(f"warning: {err}", file=sys.stderr)
    covering = [c for c in candidates
                if c[1] <= when - timedelta(minutes=10) and c[2] >= when + timedelta(minutes=10)]
    if not covering:
        raise RuntimeError("no POEORB file covers that time "
                           "(precise orbits appear ~3 weeks after acquisition)")
    url = covering[-1][0]
    out.mkdir(parents=True, exist_ok=True)
    name = url.rsplit("/", 1)[1]
    eof = out / name.removesuffix(".zip")
    if eof.exists():
        return eof
    archive = download(url, out / name)
    if archive.suffix == ".zip":
        with zipfile.ZipFile(archive) as z:
            z.extractall(out)
        archive.unlink()
    return eof


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--mission", default="S1A", choices=["S1A", "S1B", "S1C"])
    p.add_argument("--time", required=True, help="UTC, e.g. 2024-06-01T05:30:00")
    p.add_argument("--out", type=Path, default=Path("data/raw"))
    args = p.parse_args()
    try:
        eof = fetch_orbit(args.mission, datetime.fromisoformat(args.time), args.out)
    except RuntimeError as err:
        print(err, file=sys.stderr)
        return 1
    print(f"wrote {eof}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
