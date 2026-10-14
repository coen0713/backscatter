"""Apply the pre-registered scene-selection rule to the public catalogue.

docs/prereg.md, section 2: for each direction, among acquisitions between
2024-07-15 and 2024-09-15 whose footprint covers at least 95% of the region,
take the closest to 2024-08-15 that passes the snow check (validation/snow_check.py,
MeteoSwiss open data; otherwise the next closest). Tiles outside a selected scene's footprint are dropped by the
pre-registered tile exclusions. (The draft rule required full coverage; no
ascending acquisition in the window achieves that, the best covering 97.2%.)

An acquisition is one Sentinel-1 datatake: IW GRD products are cut into ~25 s
slices, so its footprint is the union of its slices' footprints. Coverage is
tested on a dense grid of points over the region (boundary and interior).

Uses only the public Copernicus Data Space catalogue (no account) and reads
metadata only, never imagery.

  python -m validation.select_scenes                       # ranked candidates
  python -m validation.select_scenes --selected build/validation/selected_scenes.json
"""

from __future__ import annotations

import argparse
import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import defaultdict
from datetime import datetime, timezone

from validation import regions, snow_check

CATALOGUE = "https://catalogue.dataspace.copernicus.eu/odata/v1/Products"
WINDOW = ("2024-07-15T00:00:00.000Z", "2024-09-15T00:00:00.000Z")
TARGET = datetime(2024, 8, 15, tzinfo=timezone.utc)
PRODUCT_TYPE = "IW_GRDH_1S"
PLATFORM = "A"  # Sentinel-1A (prereg, section 2)
MIN_COVERAGE = 0.95
# Region bounds (W, S, E, N); the snow-check station set has the same name.
REGIONS = {"bernese": (regions.LON_RANGE[0], regions.LAT_RANGE[0],  # round 1, docs/prereg.md
                       regions.LON_RANGE[1], regions.LAT_RANGE[1]),
           "graubuenden": (9.0, 46.0, 11.0, 47.0)}                  # round 2, docs/prereg_angular.md


def query(bounds: tuple[float, float, float, float] = REGIONS["bernese"],
          retries: int = 8) -> list[dict]:
    w, s, e, n = bounds
    flt = (f"Collection/Name eq 'SENTINEL-1' and "
           f"Attributes/OData.CSC.StringAttribute/any(att:att/Name eq 'productType' and "
           f"att/OData.CSC.StringAttribute/Value eq '{PRODUCT_TYPE}') and "
           f"OData.CSC.Intersects(area=geography'SRID=4326;POLYGON(({w} {s},{e} {s},{e} {n},"
           f"{w} {n},{w} {s}))') and ContentDate/Start gt {WINDOW[0]} and "
           f"ContentDate/Start lt {WINDOW[1]}")
    products, skip = [], 0
    while True:
        params = urllib.parse.urlencode({"$filter": flt, "$top": 100, "$skip": skip,
                                         "$orderby": "ContentDate/Start", "$expand": "Attributes"})
        for attempt in range(retries):
            try:
                with urllib.request.urlopen(f"{CATALOGUE}?{params}", timeout=60) as resp:
                    page = json.load(resp)["value"]
                break
            except (urllib.error.URLError, ConnectionError, TimeoutError):
                if attempt == retries - 1:
                    raise
                time.sleep(min(2.0 ** attempt, 30.0))
        products += page
        if len(page) < 100:
            return products
        skip += 100


def polygons(geojson: dict) -> list[list[tuple[float, float]]]:
    if geojson["type"] == "Polygon":
        return [[tuple(p) for p in geojson["coordinates"][0]]]
    if geojson["type"] == "MultiPolygon":
        return [[tuple(p) for p in poly[0]] for poly in geojson["coordinates"]]
    raise ValueError(f"unsupported footprint type {geojson['type']}")


def inside(pt: tuple[float, float], ring: list[tuple[float, float]]) -> bool:
    x, y = pt
    hit = False
    for (x0, y0), (x1, y1) in zip(ring, ring[1:] + ring[:1]):
        if (y0 > y) != (y1 > y) and x < x0 + (y - y0) * (x1 - x0) / (y1 - y0):
            hit = not hit
    return hit


def region_points(bounds: tuple[float, float, float, float] = REGIONS["bernese"],
                  step: float = 0.02) -> list[tuple[float, float]]:
    w, s, e, n = bounds
    xs = [w + i * step for i in range(round((e - w) / step) + 1)]
    ys = [s + j * step for j in range(round((n - s) / step) + 1)]
    return [(x, y) for x in xs for y in ys]


def acquisitions(products: list[dict],
                 bounds: tuple[float, float, float, float] = REGIONS["bernese"]) -> list[dict]:
    groups: dict[str, list[dict]] = defaultdict(list)
    for p in products:
        attrs = {a["Name"]: a.get("Value") for a in p.get("Attributes", [])}
        if attrs.get("platformSerialIdentifier") != PLATFORM:
            continue
        datatake = p["Name"].split("_")[7]  # e.g. 06BAD3
        groups[datatake].append({**p, "_attrs": attrs})
    pts = region_points(bounds)
    out = []
    for datatake, slices in groups.items():
        rings = [r for s in slices for r in polygons(s["GeoFootprint"])]
        covered = sum(any(inside(pt, r) for r in rings) for pt in pts) / len(pts)
        start = min(datetime.fromisoformat(s["ContentDate"]["Start"].replace("Z", "+00:00"))
                    for s in slices)
        a = slices[0]["_attrs"]
        out.append({"datatake": datatake, "start": start.isoformat(),
                    "direction": a.get("orbitDirection"),
                    "relative_orbit": a.get("relativeOrbitNumber"),
                    "coverage": covered,
                    "days_from_target": abs((start - TARGET).total_seconds()) / 86400.0,
                    "products": sorted(s["Name"] for s in slices),
                    "ids": [s["Id"] for s in sorted(slices, key=lambda s: s["Name"])],
                    "size_bytes": sum(int(s.get("ContentLength", 0)) for s in slices)})
    return out


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--json", help="write all ranked eligible candidates here")
    p.add_argument("--selected", help="write the selected acquisition per direction here "
                   "(input for data/fetch_s1.py --selected)")
    p.add_argument("--region", choices=sorted(REGIONS), default="bernese",
                   help="bernese = round 1 (default), graubuenden = round 2")
    args = p.parse_args()
    bounds = REGIONS[args.region]
    acqs = acquisitions(query(bounds), bounds)
    result = {}
    selected = {}
    for direction in ("ASCENDING", "DESCENDING"):
        cands = sorted((a for a in acqs if a["direction"] == direction),
                       key=lambda a: (a["coverage"] < MIN_COVERAGE, a["days_from_target"]))
        full = [a for a in cands if a["coverage"] >= MIN_COVERAGE]
        # Snow check in rank order; the first acquisition that passes is selected.
        for a in full:
            a["snow_check"] = snow_check.check(datetime.fromisoformat(a["start"]), args.region)
            if a["snow_check"]["passed"]:
                selected[direction] = a
                break
        result[direction] = full
        print(f"\n{direction}: {len(cands)} acquisitions, {len(full)} cover "
              f">= {MIN_COVERAGE:.0%} of the region")
        print("| Rank | Start (UTC) | Rel. orbit | Datatake | Days from 08-15 | Coverage | Slices | "
              "Size | Snow check |")
        print("|---|---|---|---|---|---|---|---|---|")
        for i, a in enumerate(cands[:8], start=1):
            rank = str(i) if a["coverage"] >= MIN_COVERAGE else "-"
            sc = a.get("snow_check")
            verdict = ("not needed" if sc is None and a["coverage"] >= MIN_COVERAGE else
                       "-" if sc is None else "pass" if sc["passed"] else
                       "FAIL: " + "; ".join(sc["reasons"]))
            print(f"| {rank} | {a['start'][:19]} | {a['relative_orbit']} | {a['datatake']} | "
                  f"{a['days_from_target']:.1f} | {a['coverage']:.1%} | {len(a['products'])} | "
                  f"{a['size_bytes'] / 1e9:.1f} GB | {verdict} |")
        if direction in selected:
            print("Selected: " + ", ".join(selected[direction]["products"]))
        else:
            print("No acquisition passes coverage and snow check.")
    if args.json:
        with open(args.json, "w") as f:
            json.dump(result, f, indent=1)
    if args.selected:
        with open(args.selected, "w") as f:
            json.dump(selected, f, indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
