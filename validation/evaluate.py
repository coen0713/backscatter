"""Real-data evaluation (docs/prereg.md, sections 4 and 6).

  # 1. One ellipsoidal DEM for the whole region, for SNAP's terrain correction:
  python -m validation.evaluate region-dem --out build/validation/region_dem.tif

  # 2. Run validation/snap/s1_grd_beta0_tc.xml on each scene with that DEM.

  # 3. Simulate and compare every tile of a set against the SNAP outputs:
  python -m validation.evaluate run --set calibration \
      --asc snap_asc.tif --asc-time 2024-08-14T17:15:51 \
      --desc snap_desc.tif --desc-time 2024-08-18T05:35:49 \
      --build build/release --out build/validation/calibration

`--set held-out` is refused until docs/prereg.md is registered.
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import rasterio
from rasterio.enums import Resampling

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "data"))

from fetch_orbit import fetch_orbit  # noqa: E402
from provenance import exe, header  # noqa: E402

from validation import metrics, pipeline, prepare, regions  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent

# Pre-registered thresholds (docs/prereg.md, section 4).
THRESHOLDS = {"G1": 0.60, "G2": 0.60, "G3": 0.85, "R1": 0.50, "L1": 1.0}
BARE = 60
GRASS = 30
TREES = 10


def load_snap_reference(path: Path, grid: prepare.Grid, beta0_band: str, mask_band: str,
                        layover_bit: int, shadow_bit: int) -> dict:
    """SNAP terrain-corrected GeoTIFF -> reference arrays on `grid`. Bands are
    found by description, or by 1-based index if given as a number."""
    with rasterio.open(path) as ds:
        names = list(ds.descriptions)

        def index(band: str) -> int:
            if band.isdigit():
                return int(band)
            for i, d in enumerate(names, start=1):
                if d and band.lower() in d.lower():
                    return i
            raise KeyError(f"band '{band}' not in {path.name} (descriptions: {names}); "
                           "pass a 1-based index instead")

        src = prepare.Grid(ds.transform, ds.width, ds.height, ds.crs)
        beta0 = ds.read(index(beta0_band)).astype(np.float64)
        mask = ds.read(index(mask_band)).astype(np.float64)
    # Footprint: thermal-noise removal leaves radar shadow at exactly zero (or
    # below), so "beta0 > 0" cannot define the swath. Per row, everything
    # between the first and last positive pixel is inside it; zero or negative
    # values inside are kept and clamped to the -60 dB floor.
    finite = np.isfinite(beta0)
    positive = finite & (beta0 > 0)
    footprint = np.zeros_like(positive)
    for i, row in enumerate(positive):
        idx = np.flatnonzero(row)
        if idx.size:
            footprint[i, idx[0]:idx[-1] + 1] = True
    footprint &= finite
    beta0 = np.where(footprint, np.maximum(beta0, 1e-6), prepare.NODATA)
    beta0_g = prepare.resample(beta0, src, grid, Resampling.bilinear, prepare.NODATA)
    # 0 is a legitimate mask value ("neither"), so the mask never defines validity.
    mask_g = prepare.resample(mask, src, grid, Resampling.nearest, None)
    valid = (beta0_g != prepare.NODATA) & (mask_g != prepare.NODATA)
    m = np.where(valid, mask_g, 0).astype(np.int64)
    return {"ref_beta0": np.where(valid, beta0_g, 0.0),
            "ref_layover": (m & layover_bit) != 0,
            "ref_shadow": (m & shadow_bit) != 0,
            "ref_valid": valid}


def pooled(results: list[dict]) -> dict:
    """Pre-registered pooling: masks and blocks over all kept tiles; L1 as the
    median over tiles."""
    keep = [r for r in results if not r["exclusions"]]
    mount = [r["_arrays"] for r in keep if r["_arrays"]["mountainous"]]
    out: dict = {"tiles_total": len(results), "tiles_kept": len(keep),
                 "tiles_mountainous": len(mount)}
    if mount:
        cat = {k: np.concatenate([a[k].ravel() for a in mount])
               for k in ("valid", "sim_layover", "ref_layover", "sim_shadow", "ref_shadow")}
        out["G1_layover_iou"] = metrics.iou(cat["sim_layover"], cat["ref_layover"], cat["valid"])
        out["G2_shadow_iou"] = metrics.iou(cat["sim_shadow"], cat["ref_shadow"], cat["valid"])
        # Tolerant F1 needs 2-D neighbourhoods: pool per-tile counts.
        counts = {}
        for kind in ("layover", "shadow"):
            sa = sb = pa = rb = 0
            for a in mount:
                va = a[f"sim_{kind}"] & a["valid"]
                vb = a[f"ref_{kind}"] & a["valid"]
                sa += np.count_nonzero(va)
                sb += np.count_nonzero(vb)
                pa += np.count_nonzero(va & metrics.dilate(vb))
                rb += np.count_nonzero(vb & metrics.dilate(va))
            p = pa / sa if sa else 0.0
            r = rb / sb if sb else 0.0
            counts[kind] = 2 * p * r / (p + r) if p + r else (float("nan") if not sa + sb else 0.0)
        out["G1_layover_f1"] = counts["layover"]
        out["G2_shadow_f1"] = counts["shadow"]
        out["G3_shadow_auc"] = metrics.auc_lower_is_positive(
            np.concatenate([a["g3_pos"] for a in mount]),
            np.concatenate([a["g3_neg"] for a in mount]))
    blk_sim = np.concatenate([r["_arrays"]["blk_sim"] for r in keep]) if keep else np.array([])
    blk_real = np.concatenate([r["_arrays"]["blk_real"] for r in keep]) if keep else np.array([])
    blk_cls = np.concatenate([r["_arrays"]["blk_class"] for r in keep]) if keep else np.array([])
    out["R_per_class"] = {metrics.WORLDCOVER.get(c, str(c)): v for c, v in
                          metrics.per_class_pearson(blk_sim, blk_real, blk_cls,
                                                    np.ones(blk_cls.shape, bool)).items()}
    offsets = np.array([r["L1_offset_px"] for r in keep]) if keep else np.zeros((0, 2))
    out["L1_median_offset_px"] = tuple(np.median(offsets, axis=0)) if len(offsets) else (np.nan,) * 2
    return out


def verdicts(p: dict) -> dict:
    r1 = p["R_per_class"].get(metrics.WORLDCOVER[BARE], (np.nan, 0))[0]
    l1 = p["L1_median_offset_px"]
    return {
        "G1 (layover, tolerant F1)": p.get("G1_layover_f1", np.nan) >= THRESHOLDS["G1"],
        "G2 (shadow, tolerant F1)": p.get("G2_shadow_f1", np.nan) >= THRESHOLDS["G2"],
        "G3 (shadow AUC)": p.get("G3_shadow_auc", np.nan) >= THRESHOLDS["G3"],
        "R1 (bare/sparse r)": r1 >= THRESHOLDS["R1"],
        "L1 (median offset)": bool(np.all(np.abs(l1) <= THRESHOLDS["L1"])),
    }


def run(args: argparse.Namespace) -> int:
    tiles = regions.calibration_tiles() if args.set == "calibration" else regions.held_out_tiles()
    if args.tiles:
        tiles = [t for t in tiles if t.name in args.tiles]
    for t in tiles:
        pipeline.check_allowed(t)
    render = exe(args.build, "apps", "bsar_render")
    report: dict = {}
    for direction, snap, utc in (("ascending", args.asc, args.asc_time),
                                 ("descending", args.desc, args.desc_time)):
        if snap is None:
            continue
        orbit = fetch_orbit(args.mission, datetime.fromisoformat(utc), ROOT / "data" / "raw")
        results = []
        for tile in tiles:
            sim = pipeline.simulate_tile(tile, orbit, utc, render, args.out / direction / tile.name)
            ref = load_snap_reference(snap, sim.grid, args.beta0_band, args.mask_band,
                                      args.layover_bit, args.shadow_bit)
            results.append(pipeline.compare(sim, **ref))
            print(f"{direction} {tile.name}", flush=True)
        p = pooled(results)
        report[direction] = {"pooled": p, "verdicts": verdicts(p),
                             "tiles": [{k: v for k, v in r.items() if k != "_arrays"}
                                       for r in results]}
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "results.json").write_text(json.dumps(report, indent=1, default=str))
    lines = [f"# Evaluation: {args.set} tiles\n", header("python -m validation.evaluate", args.build)]
    for direction, rep in report.items():
        p = rep["pooled"]
        lines += [f"## {direction}\n", f"Tiles: {p['tiles_kept']} kept of {p['tiles_total']} "
                  f"({p['tiles_mountainous']} mountainous).\n", "| Metric | Value |", "|---|---|"]
        for k in ("G1_layover_iou", "G1_layover_f1", "G2_shadow_iou", "G2_shadow_f1",
                  "G3_shadow_auc", "L1_median_offset_px"):
            lines.append(f"| {k} | {p.get(k, 'n/a')} |")
        for cls, (r, n) in p["R_per_class"].items():
            lines.append(f"| R, {cls} | r = {r:.3f} ({n} blocks) |")
        lines += ["", "| Criterion | Pass |", "|---|---|"]
        lines += [f"| {k} | {'yes' if v else 'NO'} |" for k, v in rep["verdicts"].items()]
        lines.append("")
    (args.out / "report.md").write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines))
    return 0


def region_dem(args: argparse.Namespace) -> int:
    m = pipeline.MARGIN_DEG
    w, s = regions.LON_RANGE[0] - m, regions.LAT_RANGE[0] - m
    e, n = regions.LON_RANGE[1] + m, regions.LAT_RANGE[1] + m
    heights, grid = prepare.ellipsoidal_dem((w, s, e, n))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with rasterio.open(args.out, "w", driver="GTiff", width=grid.width, height=grid.height,
                       count=1, dtype="float32", crs=grid.crs, transform=grid.transform,
                       nodata=prepare.NODATA, compress="deflate", tiled=True) as ds:
        ds.write(heights, 1)
    print(f"wrote {args.out} ({grid.width} x {grid.height}, ellipsoidal heights)")
    return 0


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="command", required=True)
    d = sub.add_parser("region-dem")
    d.add_argument("--out", type=Path, default=ROOT / "build" / "validation" / "region_dem.tif")
    r = sub.add_parser("run")
    r.add_argument("--set", choices=["calibration", "held-out"], required=True)
    r.add_argument("--tiles", nargs="*", help="restrict to these tile names")
    r.add_argument("--asc", type=Path)
    r.add_argument("--asc-time")
    r.add_argument("--desc", type=Path)
    r.add_argument("--desc-time")
    r.add_argument("--mission", default="S1A")
    r.add_argument("--beta0-band", default="Beta0_VV")
    r.add_argument("--mask-band", default="layover_shadow_mask")
    r.add_argument("--layover-bit", type=int, default=1)
    r.add_argument("--shadow-bit", type=int, default=2)
    r.add_argument("--build", type=Path, default=ROOT / "build" / "release")
    r.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    return region_dem(args) if args.command == "region-dem" else run(args)


if __name__ == "__main__":
    sys.exit(main())
