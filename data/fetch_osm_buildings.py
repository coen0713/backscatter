#!/usr/bin/env python3
"""Fetch OpenStreetMap building footprints for a lon/lat box via Overpass and
write them in the backscatter footprint format (crs wgs84).

Heights come from the `height` tag, else `building:levels` x 3 m, else a
default. The material comes from `building:material`/`roof:material` when it
maps onto a preset (metal, concrete), else concrete.

Example:
    python data/fetch_osm_buildings.py --bbox 8.53 47.36 8.55 47.38 \
        --out data/raw/zurich_buildings.txt
"""

from __future__ import annotations

import argparse
import json
import sys
import urllib.parse
import urllib.request
from pathlib import Path

OVERPASS = "https://overpass-api.de/api/interpreter"
LEVEL_HEIGHT = 3.0


def parse_height(tags: dict[str, str], default: float) -> float:
    for key in ("height", "building:height"):
        if key in tags:
            try:
                return float(tags[key].replace("m", "").strip())
            except ValueError:
                pass
    if "building:levels" in tags:
        try:
            return LEVEL_HEIGHT * float(tags["building:levels"])
        except ValueError:
            pass
    return default


def material(tags: dict[str, str]) -> str:
    text = " ".join(tags.get(k, "") for k in ("building:material", "roof:material")).lower()
    if "metal" in text or "steel" in text:
        return "metal"
    return "concrete"


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bbox", type=float, nargs=4, metavar=("W", "S", "E", "N"), required=True)
    p.add_argument("--default-height", type=float, default=10.0)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()

    w, s, e, n = args.bbox
    query = f"[out:json][timeout:120];way[building]({s},{w},{n},{e});out geom tags;"
    data = urllib.parse.urlencode({"data": query}).encode()
    print(f"querying Overpass for {args.bbox}")
    with urllib.request.urlopen(urllib.request.Request(OVERPASS, data=data)) as resp:
        result = json.load(resp)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with open(args.out, "w") as f:
        f.write("# OpenStreetMap building footprints (c) OpenStreetMap contributors, ODbL\n")
        f.write(f"# bbox {w} {s} {e} {n}\ncrs wgs84\n")
        for el in result.get("elements", []):
            geom = el.get("geometry")
            if el.get("type") != "way" or not geom or len(geom) < 4:
                continue
            tags = el.get("tags", {})
            f.write(f"building {parse_height(tags, args.default_height):.1f} {material(tags)}\n")
            for pt in geom[:-1]:  # OSM rings repeat the first vertex at the end
                f.write(f"{pt['lon']:.7f} {pt['lat']:.7f}\n")
            f.write("end\n")
            count += 1
    print(f"wrote {count} buildings to {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
