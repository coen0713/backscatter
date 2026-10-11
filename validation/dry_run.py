"""End-to-end dry run on calibration tiles (docs/prereg.md, section 6, step 1).

No Sentinel-1 image data is used. For each tile and pass the pipeline runs for
real (real DEM, real precise orbit, real pass geometry, WorldCover, geocoding,
UTM grid), and the metrics are exercised against synthetic stand-ins for the
reference whose answers are known:

  identity   reference = simulation           -> G1 = G2 = 1, r = 1, offset 0
  perturbed  simulation shifted by (+1, +2) px, 9-look speckle, -22 dB noise
             floor                             -> L1 recovers (+1, +2)
  cross-pass ascending vs. descending masks    -> layover on opposite slopes

Usage:
  python -m validation.dry_run --build build/release [--tiles N46.60_E007.90 ...]
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "data"))

from fetch_orbit import fetch_orbit  # noqa: E402
from provenance import exe, header  # noqa: E402

import rasterio  # noqa: E402

from validation import evaluate, metrics, pipeline, regions  # noqa: E402
from validation.quicklook import overlay, write_png  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent

# Real Sentinel-1A IW-like passes over the region (bsar_eval passes, POEORB),
# the closest to 2024-08-15 in each direction. Orbit geometry only.
DEFAULT_PASSES = {
    "ascending": "2024-08-14T17:15:51",
    "descending": "2024-08-18T05:35:49",
}
# Lauterbrunnen / Jungfrau calibration block (4, 3): the steepest terrain in the region.
DEFAULT_TILES = ["N46.60_E007.80", "N46.60_E007.90", "N46.70_E007.80", "N46.70_E007.90"]


def tile_by_name(name: str) -> regions.Tile:
    for t in regions.calibration_tiles() + regions.held_out_tiles():
        if t.name == name:
            return t
    raise ValueError(f"unknown tile {name}")


def perturbed_reference(sim: pipeline.TileProducts, rng: np.random.Generator,
                        shift: tuple[int, int] = (1, 2)) -> dict:
    beta0 = np.roll(sim.beta0, shift, axis=(0, 1))
    beta0 = beta0 * rng.gamma(9.0, 1.0 / 9.0, beta0.shape) + 10 ** (-22 / 10)  # speckle + NESZ
    return {"ref_beta0": beta0,
            "ref_layover": np.roll(sim.layover, shift, axis=(0, 1)),
            "ref_shadow": np.roll(sim.shadow, shift, axis=(0, 1)),
            "ref_valid": sim.valid.copy()}


def fmt(x) -> str:
    if isinstance(x, float):
        return "n/a" if np.isnan(x) else f"{x:.3f}"
    if isinstance(x, tuple) and len(x) == 2 and all(isinstance(v, float) for v in x):
        return f"({x[0]:+.2f}, {x[1]:+.2f})"
    return str(x)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--build", type=Path, default=ROOT / "build" / "release")
    p.add_argument("--work", type=Path, default=ROOT / "build" / "validation_dry_run")
    p.add_argument("--tiles", nargs="*", default=DEFAULT_TILES)
    p.add_argument("--mission", default="S1A")
    args = p.parse_args()

    render = exe(args.build, "apps", "bsar_render")
    tiles = [tile_by_name(n) for n in args.tiles]
    for t in tiles:
        pipeline.check_allowed(t)  # refuses held-out tiles before registration
    rng = np.random.default_rng(20261010)
    sims: dict[tuple[str, str], pipeline.TileProducts] = {}
    rows = []
    t_start = time.time()

    for direction, utc in DEFAULT_PASSES.items():
        orbit = fetch_orbit(args.mission, datetime.fromisoformat(utc), ROOT / "data" / "raw")
        for tile in tiles:
            t0 = time.time()
            work = args.work / direction / tile.name
            sim = pipeline.simulate_tile(tile, orbit, utc, render, work)
            sims[(direction, tile.name)] = sim
            write_png(work / "quicklook.png", overlay(sim.beta0, sim.layover, sim.shadow, sim.valid))
            ident = pipeline.compare(sim, sim.beta0, sim.layover, sim.shadow, sim.valid)
            pert = pipeline.compare(sim, **perturbed_reference(sim, rng))
            rows.append((direction, tile, sim, ident, pert, time.time() - t0))
            print(f"{direction} {tile.name}: {time.time() - t0:.1f} s", flush=True)

    # SNAP-format round trip: each simulation written as a SNAP-style GeoTIFF
    # (Beta0_VV + layover_shadow_mask bit flags), read back through the real
    # loader, compared, and pooled. Must reproduce the identity result.
    ok = True
    snap_results: dict[str, list[dict]] = {}
    for direction, tile, sim, _, _, _ in rows:
        fake = args.work / direction / tile.name / "fake_snap.tif"
        mask = sim.layover.astype(np.uint8) * 1 + sim.shadow.astype(np.uint8) * 2
        # A real SNAP scene is a 250 km swath, so a tile is always deep inside
        # it; embed the tile in a border of ordinary pixels to match.
        pad = 8
        beta0 = np.pad(np.where(sim.valid, sim.beta0, 0), pad, constant_values=0.1)
        mask = np.pad(mask, pad, constant_values=0)
        t = sim.grid.transform
        transform = rasterio.Affine(t.a, 0, t.c - pad * t.a, 0, t.e, t.f - pad * t.e)
        with rasterio.open(fake, "w", driver="GTiff", width=beta0.shape[1], height=beta0.shape[0],
                           count=2, dtype="float32", crs=sim.grid.crs, transform=transform,
                           nodata=0) as ds:
            ds.write(beta0.astype(np.float32), 1)
            ds.write(mask.astype(np.float32), 2)
            ds.set_band_description(1, "Beta0_VV")
            ds.set_band_description(2, "layover_shadow_mask")
        ref = evaluate.load_snap_reference(fake, sim.grid, "Beta0_VV", "layover_shadow_mask", 1, 2)
        snap_results.setdefault(direction, []).append(pipeline.compare(sim, **ref))
    snap_pooled = {d: evaluate.pooled(r) for d, r in snap_results.items()}
    for d, pp in snap_pooled.items():
        ok = ok and pp["G1_layover_f1"] == 1.0 and np.allclose(pp["L1_median_offset_px"], 0, atol=0.05)

    # Cross-pass: the two look directions should put layover on different slopes.
    cross = []
    for tile in tiles:
        a, d = sims[("ascending", tile.name)], sims[("descending", tile.name)]
        both = a.valid & d.valid
        cross.append((tile.name, metrics.iou(a.layover, d.layover, both),
                      float(a.layover[both].mean()), float(d.layover[both].mean())))

    lines = ["# Validation pipeline dry run (calibration tiles, no Sentinel-1 image data)\n",
             header("python -m validation.dry_run", args.build),
             "Real inputs: Copernicus GLO-30 (ellipsoidal via EGM2008), Sentinel-1A POEORB "
             "orbits, WorldCover 2021. Passes: " +
             ", ".join(f"{k} {v} UTC" for k, v in DEFAULT_PASSES.items()) + ". "
             "References are synthetic stand-ins with known answers (see "
             "`validation/dry_run.py`); these numbers test the pipeline, not the simulator.\n",
             "## Per tile and pass\n",
             "| Pass | Tile | Relief | Layover | Shadow | Identity: G1 / G2 / G3 / L1 | "
             "Perturbed: G1 / G2 / G3 / L1 (expect L1 = (+1, +2)) | Exclusions | Time |",
             "|---|---|---|---|---|---|---|---|---|"]
    for direction, tile, sim, ident, pert, secs in rows:
        ident_ok = (ident["G1_layover_iou"] == 1.0 or np.isnan(ident["G1_layover_iou"])) and \
            np.allclose(ident["L1_offset_px"], (0, 0), atol=0.05)
        pert_ok = np.allclose(pert["L1_offset_px"], (1, 2), atol=0.3)
        ok = ok and ident_ok and pert_ok
        lines.append(
            f"| {direction} | {tile.name} | {ident['relief_m']:.0f} m | "
            f"{sim.layover[sim.valid].mean():.1%} | {sim.shadow[sim.valid].mean():.1%} | "
            f"{fmt(ident['G1_layover_iou'])} / {fmt(ident['G2_shadow_iou'])} / "
            f"{fmt(ident['G3_shadow_auc'])} / {fmt(ident['L1_offset_px'])} | "
            f"{fmt(pert['G1_layover_iou'])} / {fmt(pert['G2_shadow_iou'])} / "
            f"{fmt(pert['G3_shadow_auc'])} / {fmt(pert['L1_offset_px'])} | "
            f"{'; '.join(ident['exclusions']) or 'none'} | {secs:.1f} s |")
    lines += ["", "Per-class correlation in the identity check (should be 1.000 wherever a "
              "class has enough pure 90 m blocks):\n", "| Pass | Tile | Classes (r, blocks) |",
              "|---|---|---|"]
    for direction, tile, _, ident, _, _ in rows:
        cls = ", ".join(f"{metrics.WORLDCOVER.get(c, c)}: {r:.3f} ({n})"
                        for c, (r, n) in sorted(ident["R_per_class"].items()))
        lines.append(f"| {direction} | {tile.name} | {cls or 'none'} |")
    lines += ["", "## Sensitivity of the mask scores to misregistration\n",
              "Simulation against itself shifted by one 30 m cell east (well within the L1 "
              "threshold of 1 pixel). Strict IoU collapses for thin masks; the tolerant F1 "
              "(1-cell tolerance) does not.\n",
              "| Pass | Tile | Layover IoU | Layover F1 | Shadow IoU | Shadow F1 |",
              "|---|---|---|---|---|---|"]
    for direction, tile, sim, _, _, _ in rows:
        moved_l = np.roll(sim.layover, 1, axis=1)
        moved_s = np.roll(sim.shadow, 1, axis=1)
        lines.append(
            f"| {direction} | {tile.name} | {fmt(metrics.iou(sim.layover, moved_l, sim.valid))} | "
            f"{fmt(metrics.tolerant_f1(sim.layover, moved_l, sim.valid))} | "
            f"{fmt(metrics.iou(sim.shadow, moved_s, sim.valid))} | "
            f"{fmt(metrics.tolerant_f1(sim.shadow, moved_s, sim.valid))} |")
    lines += ["", "## Ascending vs. descending\n",
              "Opposite look directions should put layover on opposite slopes, so the "
              "overlap of the two layover masks should be small.\n",
              "| Tile | Layover IoU (asc vs. desc) | Layover asc | Layover desc |", "|---|---|---|---|"]
    for name, iou, fa, fd in cross:
        lines.append(f"| {name} | {iou:.3f} | {fa:.1%} | {fd:.1%} |")
    lines += ["", "## SNAP-format round trip (pooled over tiles)\n",
              "Simulation written as a SNAP-style GeoTIFF and read back through "
              "`evaluate.load_snap_reference`, then pooled with `evaluate.pooled`. Expect G1/G2 "
              "= 1 and zero offset.\n",
              "| Pass | Tiles kept | G1 IoU / F1 | G2 IoU / F1 | G3 AUC | L1 median offset |",
              "|---|---|---|---|---|---|"]
    for d, pp in snap_pooled.items():
        lines.append(f"| {d} | {pp['tiles_kept']} of {pp['tiles_total']} | "
                     f"{fmt(pp['G1_layover_iou'])} / {fmt(pp['G1_layover_f1'])} | "
                     f"{fmt(pp['G2_shadow_iou'])} / {fmt(pp['G2_shadow_f1'])} | "
                     f"{fmt(pp['G3_shadow_auc'])} | {fmt(tuple(map(float, pp['L1_median_offset_px'])))} |")
    lines += ["", f"Self-consistency checks: **{'PASS' if ok else 'FAIL'}**. "
              f"Total time {time.time() - t_start:.0f} s.", ""]
    report = args.work / "report.md"
    report.write_text("\n".join(lines), encoding="utf-8")
    (args.work / "summary.json").write_text(json.dumps(
        [{"pass": d, "tile": t.name, "identity": {k: str(v) for k, v in i.items()},
          "perturbed": {k: str(v) for k, v in q.items()}} for d, t, _, i, q, _ in rows], indent=1))
    print("\n".join(lines))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
