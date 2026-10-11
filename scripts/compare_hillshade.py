#!/usr/bin/env python3
"""Compare bsar_render's ray-traced hillshade with `gdaldem hillshade` on a real DEM.

Phase 1's "done when" criterion: the ray-tracing core renders an optical
hillshade of a real DEM tile that matches GDAL's. Steps:

1. Read a crop of a Copernicus GLO-30 tile straight from AWS Open Data
   (/vsicurl/, only the needed blocks) and warp it to UTM at 30 m, so both
   tools see the same metric grid.
2. `gdaldem hillshade` (Horn's 3x3 gradient, azimuth 315, elevation 45).
3. `bsar_render --mode hillshade` on the same grid: one vertical ray per cell
   centre, shaded with the normal of the triangle it hits.
4. Compare at our pixel centres, which are the DEM cell centres, by averaging
   the four surrounding gdaldem pixels.

The two are not expected to agree exactly: Horn's operator smooths the gradient
over 3x3 samples (60 m), while the ray tracer shades individual 30 m triangle
facets. The thresholds below were fixed before the first run.

Needs the GDAL command-line tools and NumPy.
"""

from __future__ import annotations

import argparse
import shutil
import struct
import subprocess
import sys
import zlib
from pathlib import Path

import numpy as np

from provenance import ROOT, exe, header

TILE = ("https://copernicus-dem-30m.s3.amazonaws.com/Copernicus_DSM_COG_10_N46_00_E007_00_DEM/"
        "Copernicus_DSM_COG_10_N46_00_E007_00_DEM.tif")
# Lauterbrunnen valley and the Jungfrau massif: 3,000 m of relief, cliffs.
BBOX = (7.80, 46.50, 8.00, 46.62)  # W S E N, degrees
UTM_EPSG = "EPSG:32632"
RESOLUTION = 30.0

# Pre-stated pass criteria (fixed before the first comparison was run).
MIN_PEARSON_R = 0.95
MAX_ABS_BIAS = 0.02


def run(cmd: list[str]) -> None:
    print("+", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)


def read_asc(path: Path) -> tuple[np.ndarray, float | None]:
    """ESRI ASCII grid -> (array with row 0 = north, nodata value)."""
    meta: dict[str, float] = {}
    skip = 0
    with open(path) as f:
        for line in f:
            key, *rest = line.split()
            if not rest or key[0].isdigit() or key[0] in "-.":
                break
            meta[key.lower()] = float(rest[0])
            skip += 1
    data = np.loadtxt(path, skiprows=skip, dtype=np.float64)
    return data.reshape(int(meta["nrows"]), int(meta["ncols"])), meta.get("nodata_value")


def write_png(path: Path, gray: np.ndarray) -> None:
    img = np.clip(gray, 0, 255).astype(np.uint8)
    raw = b"".join(b"\x00" + row.tobytes() for row in img)

    def chunk(kind: bytes, body: bytes) -> bytes:
        return struct.pack(">I", len(body)) + kind + body + struct.pack(">I", zlib.crc32(kind + body))

    ihdr = struct.pack(">IIBBBBB", img.shape[1], img.shape[0], 8, 0, 0, 0, 0)
    path.write_bytes(b"\x89PNG\r\n\x1a\n" + chunk(b"IHDR", ihdr) +
                     chunk(b"IDAT", zlib.compress(raw, 9)) + chunk(b"IEND", b""))


def correlation(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.corrcoef(a.ravel(), b.ravel())[0, 1])


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--build", type=Path, default=ROOT / "build" / "release")
    p.add_argument("--work", type=Path, default=ROOT / "build" / "hillshade_compare")
    p.add_argument("--dem", help="use this GDAL-readable DEM instead of the Copernicus crop")
    p.add_argument("--strict", action="store_true", help="exit non-zero if the criteria fail")
    args = p.parse_args()

    for tool in ("gdalwarp", "gdal_translate", "gdaldem"):
        if shutil.which(tool) is None:
            print(f"{tool} not found: install the GDAL command-line tools", file=sys.stderr)
            return 2
    work = args.work
    work.mkdir(parents=True, exist_ok=True)
    render = exe(args.build, "apps", "bsar_render")

    # 1. Same metric grid for both tools.
    dem_tif = work / "dem_utm.tif"
    src = args.dem or f"/vsicurl/{TILE}"
    w, s, e, n = BBOX
    run(["gdalwarp", "-q", "-overwrite", "-t_srs", UTM_EPSG, "-tr", str(RESOLUTION), str(RESOLUTION),
         "-r", "bilinear", "-te_srs", "EPSG:4326", "-te", str(w), str(s), str(e), str(n),
         "-ot", "Float32", src, str(dem_tif)])
    dem_asc = work / "dem_utm.asc"
    run(["gdal_translate", "-q", "-of", "AAIGrid", str(dem_tif), str(dem_asc)])

    # 2. GDAL's hillshade.
    gdal_tif = work / "gdaldem_hillshade.tif"
    run(["gdaldem", "hillshade", "-q", "-az", "315", "-alt", "45", "-z", "1", "-s", "1",
         str(dem_tif), str(gdal_tif)])
    gdal_asc = work / "gdaldem_hillshade.asc"
    run(["gdal_translate", "-q", "-of", "AAIGrid", str(gdal_tif), str(gdal_asc)])

    # 3. Ours, on the same grid.
    prefix = work / "bsar"
    run([str(render), "--scene", f"dem:{dem_asc}", "--mode", "hillshade", "--pixel", str(RESOLUTION),
         "--sun-azimuth", "315", "--sun-elevation", "45", "--out", str(prefix)])
    ours = np.load(f"{prefix}_hillshade.npy").astype(np.float64)

    # 4. Align: gdaldem is per DEM sample (0 = no data, else ~1 + 254 cos);
    #    our pixels sit at the centres of the cells between samples.
    g, _ = read_asc(gdal_asc)
    dem, dem_nodata = read_asc(dem_asc)
    if ours.shape != (g.shape[0] - 1, g.shape[1] - 1):
        print(f"grid mismatch: ours {ours.shape}, gdaldem {g.shape}", file=sys.stderr)
        return 1
    valid_g = g > 0
    shade_g = np.where(valid_g, (g - 1.0) / 254.0, 0.0)
    valid = valid_g[:-1, :-1] & valid_g[1:, :-1] & valid_g[:-1, 1:] & valid_g[1:, 1:]
    if dem_nodata is not None:
        dv = dem != dem_nodata
        valid &= dv[:-1, :-1] & dv[1:, :-1] & dv[:-1, 1:] & dv[1:, 1:]
    ref = 0.25 * (shade_g[:-1, :-1] + shade_g[1:, :-1] + shade_g[:-1, 1:] + shade_g[1:, 1:])

    a = ours[valid]
    b = ref[valid]
    diff = a - b
    r = correlation(a, b)
    bias = float(diff.mean())
    mae = float(np.abs(diff).mean())
    p95 = float(np.percentile(np.abs(diff), 95))

    # Alignment sanity check: the correlation must peak at zero shift.
    shifts = {}
    core = (slice(2, -2), slice(2, -2))
    for dy in (-1, 0, 1):
        for dx in (-1, 0, 1):
            o = ours[2 + dy:ours.shape[0] - 2 + dy, 2 + dx:ours.shape[1] - 2 + dx]
            m = valid[core] & valid[2 + dy:valid.shape[0] - 2 + dy, 2 + dx:valid.shape[1] - 2 + dx]
            shifts[(dy, dx)] = correlation(o[m], ref[core][m])
    best_shift = max(shifts, key=shifts.get)

    passed = r >= MIN_PEARSON_R and abs(bias) <= MAX_ABS_BIAS and best_shift == (0, 0)

    # Side-by-side image: gdaldem | bsar | difference (mid-grey = 0, +-0.25 full scale).
    sep = np.full((ours.shape[0], 6), 255.0)
    d_img = np.where(valid, 128.0 + 512.0 * (ours - ref), 0.0)
    panel = np.hstack([np.where(valid, ref * 255.0, 0.0), sep, ours * 255.0, sep, d_img])
    write_png(work / "hillshade_compare.png", panel)

    body = [
        "# Ray-traced hillshade vs. `gdaldem hillshade`\n",
        header("scripts/compare_hillshade.py", args.build),
        f"DEM: Copernicus GLO-30, lon {w}..{e}, lat {s}..{n} (Lauterbrunnen / Jungfrau), "
        f"warped to {UTM_EPSG} at {RESOLUTION:.0f} m: {dem.shape[1]} x {dem.shape[0]} samples, "
        f"elevation {dem[dem != dem_nodata].min():.0f}..{dem[dem != dem_nodata].max():.0f} m. "
        "Sun azimuth 315 deg, elevation 45 deg, no cast shadows (matching gdaldem).\n",
        f"Compared at {valid.sum()} cell centres (gdaldem averaged over the 4 surrounding "
        "samples; edges and no-data excluded).\n",
        "| Metric | Value | Criterion (fixed before the first run) | Result |",
        "|---|---|---|---|",
        f"| Pearson r | {r:.4f} | >= {MIN_PEARSON_R} | {'pass' if r >= MIN_PEARSON_R else 'FAIL'} |",
        f"| Mean difference (bsar - gdaldem) | {bias:+.4f} | abs <= {MAX_ABS_BIAS} | "
        f"{'pass' if abs(bias) <= MAX_ABS_BIAS else 'FAIL'} |",
        f"| Mean absolute difference | {mae:.4f} | reported | |",
        f"| 95th percentile absolute difference | {p95:.4f} | reported | |",
        f"| Correlation peak at zero shift | best shift {best_shift} | (0, 0) | "
        f"{'pass' if best_shift == (0, 0) else 'FAIL'} |",
        "",
        "Correlation for +-1 cell shifts (rows: dy, columns: dx):\n",
        "| dy \\ dx | -1 | 0 | +1 |", "|---|---|---|---|",
    ]
    for dy in (-1, 0, 1):
        body.append(f"| {dy:+d} | " + " | ".join(f"{shifts[(dy, dx)]:.4f}" for dx in (-1, 0, 1)) + " |")
    body += ["", f"Overall: **{'PASS' if passed else 'FAIL'}**", "",
             "Values are in [0, 1] (cosine of the angle between surface normal and sun). "
             "Expected differences: Horn's 3x3 operator smooths the gradient over 60 m, the ray "
             "tracer shades single 30 m triangle facets, so sharp ridges and gullies differ most.",
             "", "![gdaldem | bsar | difference](images/hillshade_compare.png)", ""]
    (work / "report.md").write_text("\n".join(body))
    print("\n".join(body))
    return 1 if (args.strict and not passed) else 0


if __name__ == "__main__":
    sys.exit(main())
