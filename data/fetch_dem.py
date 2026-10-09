#!/usr/bin/env python3
"""Download a Copernicus GLO-30 DEM tile and convert it for bsar_render.

The tiles are public on the AWS Open Data registry (no account needed):
https://registry.opendata.aws/copernicus-dem/

Copernicus heights are relative to the EGM2008 geoid. Sentinel-1 orbits are
ellipsoidal (WGS84), so pass --ellipsoidal to convert with GDAL (needs the
PROJ geoid grid; `projsync --file us_nga_egm08_25` if it is missing). Without
it, absolute slant ranges are off by the local geoid height (tens of metres),
which matters when comparing against real acquisitions.

Example:
    python data/fetch_dem.py --lat 46 --lon 7 --bbox 7.85 46.45 8.05 46.60 \
        --ellipsoidal --out data/raw/eiger.asc
    bsar_render --scene dem:data/raw/eiger.asc --geographic --mode hillshade
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path

BASE = "https://copernicus-dem-30m.s3.amazonaws.com"


def tile_name(lat: int, lon: int) -> str:
    ns = "N" if lat >= 0 else "S"
    ew = "E" if lon >= 0 else "W"
    return f"Copernicus_DSM_COG_10_{ns}{abs(lat):02d}_00_{ew}{abs(lon):03d}_00_DEM"


def download(url: str, dest: Path) -> None:
    if dest.exists():
        print(f"cached: {dest}")
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"downloading {url}")
    tmp = dest.with_suffix(dest.suffix + ".part")
    with urllib.request.urlopen(url) as resp, open(tmp, "wb") as f:
        shutil.copyfileobj(resp, f)
    tmp.rename(dest)


def run(cmd: list[str]) -> None:
    print("+", " ".join(cmd))
    subprocess.run(cmd, check=True)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--lat", type=int, required=True, help="tile south edge, integer degrees")
    p.add_argument("--lon", type=int, required=True, help="tile west edge, integer degrees")
    p.add_argument("--bbox", type=float, nargs=4, metavar=("W", "S", "E", "N"),
                   help="crop to this lon/lat box")
    p.add_argument("--ellipsoidal", action="store_true",
                   help="convert EGM2008 heights to WGS84 ellipsoidal heights")
    p.add_argument("--out", type=Path, required=True, help="output .asc (or .tif) path")
    p.add_argument("--cache", type=Path, default=Path("data/cache"))
    args = p.parse_args()

    if shutil.which("gdal_translate") is None:
        print("GDAL command-line tools are required (gdal_translate, gdalwarp)", file=sys.stderr)
        return 1

    name = tile_name(args.lat, args.lon)
    tif = args.cache / f"{name}.tif"
    download(f"{BASE}/{name}/{name}.tif", tif)

    src = tif
    if args.ellipsoidal:
        ell = args.cache / f"{name}_ellipsoidal.tif"
        if not ell.exists():
            run(["gdalwarp", "-s_srs", "EPSG:4326+3855", "-t_srs", "EPSG:4979", str(tif), str(ell)])
        src = ell

    args.out.parent.mkdir(parents=True, exist_ok=True)
    cmd = ["gdal_translate"]
    if args.bbox:
        w, s, e, n = args.bbox
        cmd += ["-projwin", str(w), str(n), str(e), str(s)]
    if args.out.suffix == ".asc":
        cmd += ["-of", "AAIGrid", "-co", "FORCE_CELLSIZE=FALSE"]
    run(cmd + [str(src), str(args.out)])
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
