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
import shutil
import sys
import urllib.request
import zipfile
from datetime import datetime, timedelta
from pathlib import Path

BASE = "https://step.esa.int/auxdata/orbits/Sentinel-1/POEORB"
PATTERN = re.compile(r"(S1[ABC]_OPER_AUX_POEORB_OPOD_\d{8}T\d{6}_V(\d{8}T\d{6})_(\d{8}T\d{6})\.EOF(?:\.zip)?)")


def listing(mission: str, month: datetime) -> list[tuple[str, datetime, datetime]]:
    url = f"{BASE}/{mission}/{month:%Y}/{month:%m}/"
    with urllib.request.urlopen(url) as resp:
        html = resp.read().decode()
    out = []
    for name, start, stop in sorted(set(PATTERN.findall(html))):
        out.append((url + name, datetime.strptime(start, "%Y%m%dT%H%M%S"),
                    datetime.strptime(stop, "%Y%m%dT%H%M%S")))
    return out


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--mission", default="S1A", choices=["S1A", "S1B", "S1C"])
    p.add_argument("--time", required=True, help="UTC, e.g. 2024-06-01T05:30:00")
    p.add_argument("--out", type=Path, default=Path("data/raw"))
    args = p.parse_args()

    t = datetime.fromisoformat(args.time)
    candidates = []
    for month in {t.replace(day=1), (t + timedelta(days=2)).replace(day=1)}:
        try:
            candidates += listing(args.mission, month)
        except OSError as err:
            print(f"warning: {err}", file=sys.stderr)
    covering = [c for c in candidates if c[1] <= t - timedelta(minutes=10) and c[2] >= t + timedelta(minutes=10)]
    if not covering:
        print("no POEORB file covers that time (precise orbits appear ~3 weeks after acquisition)",
              file=sys.stderr)
        return 1
    url = covering[-1][0]
    args.out.mkdir(parents=True, exist_ok=True)
    dest = args.out / url.rsplit("/", 1)[1]
    print(f"downloading {url}")
    with urllib.request.urlopen(url) as resp, open(dest, "wb") as f:
        shutil.copyfileobj(resp, f)
    if dest.suffix == ".zip":
        with zipfile.ZipFile(dest) as z:
            z.extractall(args.out)
        dest.unlink()
        dest = args.out / dest.stem
    print(f"wrote {dest}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
