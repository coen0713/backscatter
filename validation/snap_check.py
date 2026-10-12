"""Pre-registration format check of the SNAP reference (docs/prereg.md, section 6, step 1).

Runs the production SNAP graph (validation/snap/s1_grd_beta0_tc.xml, or the
two-slice variant when two products are given) with one node added, a Subset
after terrain correction, so that one calibration block only is computed and
written; then reports what the evaluation needs to
know about SNAP's output: band names, mask values, CRS and grid. The mask
encoding is checked by physics alone, without the simulator, so no
pre-registered metric is computed on real data: layover lies on slopes that
rise steeply along the look direction, shadow on slopes that fall steeply
away from the sensor.

  python -m validation.snap_check --direction ascending --scene data/raw/S1A_..._06B719_AF33.zip
  python -m validation.snap_check --direction descending --scene data/raw/..._2F40.zip data/raw/..._8ECE.zip
"""

from __future__ import annotations

import argparse
import math
import os
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.warp import transform_bounds

from validation import prepare, regions

ROOT = Path(__file__).resolve().parent.parent
GRAPHS = {1: ROOT / "validation" / "snap" / "s1_grd_beta0_tc.xml",
          2: ROOT / "validation" / "snap" / "s1_grd_beta0_tc_2slices.xml"}
SUBSET = """  <node id="Subset">
    <operator>Subset</operator>
    <sources>
      <sourceProduct refid="Terrain-Correction"/>
    </sources>
    <parameters>
      <geoRegion>${region}</geoRegion>
      <copyMetadata>false</copyMetadata>
    </parameters>
  </node>
"""
GPT = Path(os.environ.get("SNAP_GPT", r"C:\Program Files\esa-snap\bin\gpt.exe"))
# Calibration block containing the development area (Lauterbrunnen): steep
# relief on both sides of the valley, so layover and shadow both occur.
BLOCK = (7.8, 46.6, 8.0, 46.8)
INSET = 0.01  # degrees; keeps the UTM bounding box of the subset inside the block
# Sentinel-1 is right-looking; ground-track heading over 46.5 N (from the orbit).
HEADING_DEG = {"ascending": -12.4, "descending": -167.6}


def check_block() -> tuple[float, float, float, float]:
    cal, _ = regions.split()
    if BLOCK not in [regions.block_bounds(b) for b in cal]:
        raise SystemExit(f"{BLOCK} is not a calibration block; refusing to process it")
    w, s, e, n = BLOCK
    return (w + INSET, s + INSET, e - INSET, n - INSET)


def check_graph(scenes: int, out_dir: Path) -> Path:
    """The production graph with a Subset between terrain correction and Write."""
    xml = GRAPHS[scenes].read_text(encoding="utf-8")
    write = re.compile(r'(  <node id="Write">\s*<operator>Write</operator>\s*<sources>\s*'
                       r'<sourceProduct refid=")Terrain-Correction("/>)')
    xml, count = write.subn(lambda m: SUBSET + m[1] + "Subset" + m[2], xml)
    if count != 1:
        raise SystemExit(f"unexpected Write node in {GRAPHS[scenes].name}")
    path = out_dir / f"check_{GRAPHS[scenes].name}"
    path.write_text(xml, encoding="utf-8")
    return path


def run_gpt(scenes: list[Path], dem: Path, out: Path, region: tuple[float, ...]) -> None:
    w, s, e, n = region
    wkt = f"POLYGON(({w} {s},{e} {s},{e} {n},{w} {n},{w} {s}))"
    inputs = ([f"-Pinput={scenes[0]}"] if len(scenes) == 1 else
              [f"-Pinput{i}={s}" for i, s in enumerate(scenes, start=1)])
    cmd = [str(GPT), str(check_graph(len(scenes), out.parent)), *inputs, f"-Pdem={dem}",
           f"-Poutput={out}", f"-Pregion={wkt}"]
    print(" ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)


def slope_along_look(dem_path: Path, grid: prepare.Grid, heading_deg: float) -> np.ndarray:
    """dh/dx along the horizontal look direction on the output grid (rise per metre)."""
    heights, src, _ = prepare.read_raster(dem_path)
    h = prepare.resample(heights, src, grid, Resampling.bilinear, prepare.NODATA)
    h = np.where(h == prepare.NODATA, np.nan, h)
    dy, dx = np.gradient(h, abs(grid.transform.e), grid.transform.a)
    dy = -dy  # rows run south
    look = math.radians(heading_deg + 90.0)  # azimuth, clockwise from north
    return dx * math.sin(look) + dy * math.cos(look)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--scene", type=Path, nargs="+", required=True,
                   help="one GRD product, or the two slices of one acquisition in order")
    p.add_argument("--direction", choices=sorted(HEADING_DEG), required=True)
    p.add_argument("--dem", type=Path, default=ROOT / "build" / "validation" / "region_dem.tif")
    p.add_argument("--out", type=Path, default=ROOT / "build" / "validation" / "snap_check")
    args = p.parse_args()
    region = check_block()
    args.out.mkdir(parents=True, exist_ok=True)
    tif = args.out / f"snap_{args.direction}_calblock.tif"
    if not tif.exists():
        if len(args.scene) > 2:
            raise SystemExit("one or two GRD products only")
        run_gpt([s.resolve() for s in args.scene], args.dem.resolve(), tif, region)

    with rasterio.open(tif) as ds:
        print(f"\n{tif.name}: {ds.width} x {ds.height}, {ds.count} bands, CRS {ds.crs}, "
              f"pixel {ds.transform.a:.3f} x {ds.transform.e:.3f} m, nodata {ds.nodata}")
        print(f"  descriptions: {list(ds.descriptions)}  dtypes: {list(ds.dtypes)}")
        tags = {i: ds.tags(i) for i in range(1, ds.count + 1)}
        for i, t in tags.items():
            if t:
                print(f"  band {i} tags: {t}")
        lonlat = transform_bounds(ds.crs, "EPSG:4326", *ds.bounds)
        inside = (lonlat[0] >= BLOCK[0] and lonlat[1] >= BLOCK[1]
                  and lonlat[2] <= BLOCK[2] and lonlat[3] <= BLOCK[3])
        print(f"  extent (lon/lat) {tuple(round(v, 4) for v in lonlat)}: "
              f"{'inside' if inside else 'OUTSIDE'} calibration block {BLOCK}")
        grid = prepare.Grid(ds.transform, ds.width, ds.height, ds.crs)
        bands = [ds.read(i) for i in range(1, ds.count + 1)]

    # The mask is the band holding a few small integer values; beta0 is continuous.
    for i, b in enumerate(bands, start=1):
        vals = np.unique(b[np.isfinite(b)])
        kind = "categorical" if vals.size <= 16 else "continuous"
        print(f"  band {i}: {kind}, " + (f"values {vals.tolist()}" if vals.size <= 16 else
              f"range {np.nanmin(b):.4g} .. {np.nanmax(b):.4g}, "
              f"{np.mean(b > 0):.1%} positive"))
    masks = [b for b in bands if np.unique(b[np.isfinite(b)]).size <= 16]
    if len(masks) != 1:
        print("Could not identify a single mask band.")
        return 1 if inside else 2
    mask = masks[0]
    g = slope_along_look(args.dem, grid, HEADING_DEG[args.direction])
    print("\n  Mask value vs. terrain slope along the look direction "
          "(layover: steep rise, shadow: steep fall):")
    print("  | Value | Pixels | Share | Median slope along look (deg) |")
    print("  |---|---|---|---|")
    for v in np.unique(mask[np.isfinite(mask)]):
        sel = (mask == v) & np.isfinite(g)
        med = math.degrees(math.atan(float(np.median(g[sel])))) if sel.any() else float("nan")
        print(f"  | {int(v)} | {int(sel.sum())} | {sel.mean():.1%} | {med:+.1f} |")
    return 0 if inside else 2


if __name__ == "__main__":
    sys.exit(main())
