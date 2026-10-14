"""EXPLORATORY (not pre-registered): per-class correlation with the
local-incidence trend removed, and each class's angular response.

The pre-registered R metrics (docs/prereg.md, section 4) correlate simulated
and real 90 m block means pooled over mountainous tiles, so most of their
variance is terrain. This asks what is left once local incidence is accounted
for, and whether the real angular response of each class has the shape the
model gives it (volume scattering would make it flatter).

Reads the held-out renders already made by `validation.evaluate run --set
held-out` (it never renders: a missing or stale render is an error) and the
SNAP references; same tiles (exclusions applied), same class-pure blocks as
the R metrics. Local incidence per 30 m cell is computed from the DEM normal
and the look vector at the tile centre (zero-Doppler on the POEORB orbit).

  python -m validation.explore_incidence --heldout build/validation/heldout \\
      --asc build/validation/snap_ascending.tif --asc-time 2024-08-07T17:24:01 \\
      --desc build/validation/snap_descending.tif --desc-time 2024-08-06T05:35:49 \\
      --out build/validation/explore_incidence
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import xml.etree.ElementTree as ET
from datetime import datetime
from pathlib import Path

import numpy as np
from rasterio.enums import Resampling
from rasterio.warp import transform as warp_transform

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "data"))

from fetch_orbit import fetch_orbit  # noqa: E402

from validation import evaluate, metrics, pipeline, prepare, regions  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
A, F = 6378137.0, 1 / 298.257223563  # WGS84
E2 = F * (2 - F)
BIN_DEG = 2.0
MIN_PER_BIN = 30
CLASSES = (10, 30, 40, 50, 60, 80, 100)


# ------------------------------------------------------------ orbit geometry

def read_eof(path: Path) -> tuple[np.ndarray, np.ndarray, datetime]:
    """State vectors: times (s from the first), (n, 6) ECEF position+velocity."""
    root = ET.parse(path).getroot()
    times, states = [], []
    for osv in root.iter("OSV"):
        t = datetime.fromisoformat(osv.findtext("UTC").removeprefix("UTC="))
        times.append(t)
        states.append([float(osv.findtext(k)) for k in ("X", "Y", "Z", "VX", "VY", "VZ")])
    t0 = times[0]
    return np.array([(t - t0).total_seconds() for t in times]), np.array(states), t0


def state_at(ts: np.ndarray, sv: np.ndarray, t: float) -> tuple[np.ndarray, np.ndarray]:
    """Cubic Hermite interpolation of position and velocity."""
    i = int(np.clip(np.searchsorted(ts, t) - 1, 0, len(ts) - 2))
    h = ts[i + 1] - ts[i]
    s = (t - ts[i]) / h
    p0, v0, p1, v1 = sv[i, :3], sv[i, 3:], sv[i + 1, :3], sv[i + 1, 3:]
    h00, h10, h01, h11 = 2*s**3 - 3*s**2 + 1, s**3 - 2*s**2 + s, -2*s**3 + 3*s**2, s**3 - s**2
    pos = h00 * p0 + h10 * h * v0 + h01 * p1 + h11 * h * v1
    d00, d10, d01, d11 = 6*s**2 - 6*s, 3*s**2 - 4*s + 1, -6*s**2 + 6*s, 3*s**2 - 2*s
    vel = (d00 * p0 + d01 * p1) / h + d10 * v0 + d11 * v1
    return pos, vel


def ecef(lat: float, lon: float, h: float) -> np.ndarray:
    la, lo = math.radians(lat), math.radians(lon)
    n = A / math.sqrt(1 - E2 * math.sin(la) ** 2)
    return np.array([(n + h) * math.cos(la) * math.cos(lo), (n + h) * math.cos(la) * math.sin(lo),
                     (n * (1 - E2) + h) * math.sin(la)])


def look_geometry(ts, sv, t_guess: float, lat: float, lon: float, h: float) -> tuple[float, float]:
    """(ellipsoid incidence, look azimuth clockwise from true north), degrees,
    at the zero-Doppler time of the point."""
    p = ecef(lat, lon, h)
    t = t_guess
    for _ in range(30):
        s, v = state_at(ts, sv, t)
        f = np.dot(v, p - s)
        dt = f / np.dot(v, v)  # d f/dt ~ -|v|^2
        t += dt
        if abs(dt) < 1e-7:
            break
    s, _ = state_at(ts, sv, t)
    look = (p - s) / np.linalg.norm(p - s)            # sensor -> ground
    la, lo = math.radians(lat), math.radians(lon)
    east = np.array([-math.sin(lo), math.cos(lo), 0.0])
    north = np.array([-math.sin(la) * math.cos(lo), -math.sin(la) * math.sin(lo), math.cos(la)])
    up = np.array([math.cos(la) * math.cos(lo), math.cos(la) * math.sin(lo), math.sin(la)])
    to_sensor = -look
    inc = math.degrees(math.acos(np.dot(to_sensor, up)))
    az = math.degrees(math.atan2(np.dot(look, east), np.dot(look, north)))
    return inc, az


def local_incidence(dem: np.ndarray, grid: prepare.Grid, inc_deg: float, az_true_deg: float,
                    lat: float, lon: float) -> np.ndarray:
    """Local incidence (deg) per cell from the DEM normal on the UTM grid."""
    # Grid convergence: azimuth of grid north measured from true north.
    x, y = warp_transform("EPSG:4326", grid.crs, [lon, lon], [lat, lat + 0.01])
    gamma = math.degrees(math.atan2(x[1] - x[0], y[1] - y[0]))
    az = math.radians(az_true_deg + gamma)            # look azimuth in grid coordinates
    inc = math.radians(inc_deg)
    s = np.array([-math.sin(az) * math.sin(inc), -math.cos(az) * math.sin(inc), math.cos(inc)])
    h = np.where(dem == prepare.NODATA, np.nan, dem)
    dzdy_rows, dzdx = np.gradient(h, abs(grid.transform.e), grid.transform.a)
    dzdy = -dzdy_rows                                 # rows run south
    norm = np.sqrt(dzdx ** 2 + dzdy ** 2 + 1.0)
    cos_loc = (-dzdx * s[0] - dzdy * s[1] + s[2]) / norm
    return np.degrees(np.arccos(np.clip(cos_loc, -1.0, 1.0)))


# ------------------------------------------------------------ data

def cached_render(work: Path, prep: pipeline.PreparedTile) -> pipeline.TileProducts:
    """The held-out render as evaluate made it; never renders."""
    if not (work / "sim.done").exists():
        raise SystemExit(f"no completed render in {work}; run the held-out evaluation first")
    out = {}
    for product, resampling in (("beta0", Resampling.bilinear), ("layover", Resampling.nearest),
                                ("shadow", Resampling.nearest), ("valid", Resampling.nearest)):
        values, src_grid, nodata = prepare.read_raster(work / f"sim_geo_{product}.asc")
        out[product] = prepare.resample(values, src_grid, prep.grid, resampling, nodata)
    valid = (out["valid"] > 0.5) & (out["beta0"] != prepare.NODATA) & (prep.dem != prepare.NODATA)
    return pipeline.TileProducts(prep.tile, prep.grid, np.where(valid, out["beta0"], 0.0),
                                 (out["layover"] > 0.5) & valid, (out["shadow"] > 0.5) & valid,
                                 valid, prep.dem, prep.worldcover)


def tile_blocks(sim, ref, theta) -> dict[str, np.ndarray]:
    """Class-pure 90 m blocks exactly as pipeline.compare builds them for R1-R3,
    plus the block-mean local incidence."""
    valid = sim.valid & ref["ref_valid"]
    real_lin, ok = metrics.block_reduce(ref["ref_beta0"], valid, pipeline.BLOCK)
    sim_lin, _ = metrics.block_reduce(sim.beta0, valid, pipeline.BLOCK)
    th, _ = metrics.block_reduce(np.nan_to_num(theta, nan=0.0), valid, pipeline.BLOCK)
    classes = np.zeros(real_lin.shape, dtype=np.int32)
    for c in np.unique(sim.worldcover):
        if c in (0, pipeline.SNOW_ICE):
            continue
        frac = metrics.block_fraction(sim.worldcover == c, valid, pipeline.BLOCK)
        classes[frac * pipeline.BLOCK ** 2 >= pipeline.CLASS_PURITY - 1e-9] = c
    sel = ok & (classes > 0)
    return {"sim": metrics.to_db(sim_lin)[sel], "real": metrics.to_db(real_lin)[sel],
            "cls": classes[sel], "theta": th[sel]}


# ------------------------------------------------------------ analysis

def trend(theta: np.ndarray, y: np.ndarray, edges: np.ndarray) -> np.ndarray:
    """Binned-mean curve of y against theta, linearly interpolated (NaN outside
    bins with fewer than MIN_PER_BIN blocks)."""
    idx = np.digitize(theta, edges) - 1
    centres, means = [], []
    for b in range(len(edges) - 1):
        m = idx == b
        if m.sum() >= MIN_PER_BIN:
            centres.append(0.5 * (edges[b] + edges[b + 1]))
            means.append(y[m].mean())
    if len(centres) < 2:
        return np.full_like(y, np.nan)
    out = np.interp(theta, centres, means)
    out[(theta < centres[0]) | (theta > centres[-1])] = np.nan
    return out


def analyse(d: dict[str, np.ndarray]) -> dict:
    edges = np.arange(0.0, 90.0 + BIN_DEG, BIN_DEG)
    res = {}
    for c in CLASSES:
        m = d["cls"] == c
        if m.sum() < 50:
            continue
        th, s, r = d["theta"][m], d["sim"][m], d["real"][m]
        ts, tr = trend(th, s, edges), trend(th, r, edges)
        ok = np.isfinite(ts) & np.isfinite(tr)
        if ok.sum() < 50:
            continue
        rs, rr = s[ok] - ts[ok], r[ok] - tr[ok]
        r2 = lambda y, t: 1.0 - np.var(y - t) / np.var(y)  # noqa: E731
        # Angular response: linear slope of dB against local incidence, 20-60 deg.
        band = ok & (th >= 20) & (th <= 60)
        slope = lambda y: float(np.polyfit(th[band], y[band], 1)[0] * 10) if band.sum() > 50 else np.nan  # noqa: E731
        res[c] = {"blocks": int(m.sum()), "blocks_used": int(ok.sum()),
                  "r_raw": float(np.corrcoef(s, r)[0, 1]),
                  "r_detrended": float(np.corrcoef(rs, rr)[0, 1]),
                  "r2_incidence_real": float(r2(r[ok], tr[ok])),
                  "r2_incidence_sim": float(r2(s[ok], ts[ok])),
                  "slope_real_db_per_10deg": slope(r), "slope_sim_db_per_10deg": slope(s),
                  "theta_p5_p95": (float(np.percentile(th, 5)), float(np.percentile(th, 95))),
                  "curve": curve_table(th, s, r)}
    return res


def curve_table(th, s, r, step=10.0) -> list[tuple]:
    rows = []
    for lo in np.arange(0.0, 90.0, step):
        m = (th >= lo) & (th < lo + step)
        if m.sum() >= MIN_PER_BIN:
            rows.append((float(lo), int(m.sum()), float(np.median(r[m])), float(np.median(s[m]))))
    return rows


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--heldout", type=Path, required=True, help="output of the held-out evaluation")
    p.add_argument("--asc", type=Path, required=True)
    p.add_argument("--asc-time", required=True)
    p.add_argument("--desc", type=Path, required=True)
    p.add_argument("--desc-time", required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--tiles", nargs="*", help="restrict to these kept tiles (for a quick check)")
    args = p.parse_args()
    results = json.loads((args.heldout / "results.json").read_text())
    kept = sorted({t["tile"] for t in results["ascending"]["tiles"] if not t["exclusions"]})
    if args.tiles:
        kept = [t for t in kept if t in args.tiles]
    held = {t.name: t for t in regions.held_out_tiles()}
    report = {}
    for direction, snap, utc in (("ascending", args.asc, args.asc_time),
                                 ("descending", args.desc, args.desc_time)):
        when = datetime.fromisoformat(utc)
        ts, sv, t0 = read_eof(fetch_orbit("S1A", when, ROOT / "data" / "raw"))
        t_guess = (when - t0).total_seconds()
        parts, geom = [], {}
        for name in kept:
            tile = held[name]
            work = args.heldout / direction / name
            prep = pipeline.prepare_tile(tile, args.out / "prep" / direction / name)
            sim = cached_render(work, prep)
            ref = evaluate.load_snap_reference(snap, prep.grid, "1", "2", 1, 2)
            w, s, e, n = tile.bounds
            lat, lon = 0.5 * (s + n), 0.5 * (w + e)
            h = float(np.nanmedian(np.where(prep.dem == prepare.NODATA, np.nan, prep.dem)))
            inc, az = look_geometry(ts, sv, t_guess, lat, lon, h)
            geom[name] = (round(inc, 2), round(az, 2))
            theta = local_incidence(prep.dem, prep.grid, inc, az, lat, lon)
            parts.append(tile_blocks(sim, ref, theta))
            print(f"{direction} {name}: incidence {inc:.2f}, look azimuth {az:.1f}", flush=True)
        d = {k: np.concatenate([p_[k] for p_ in parts]) for k in parts[0]}
        report[direction] = {"classes": analyse(d), "tile_geometry": geom}
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "results.json").write_text(json.dumps(report, indent=1))
    lines = render_report(report, len(kept))
    (args.out / "report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    return 0


def render_report(report: dict, n_tiles: int) -> list[str]:
    name = metrics.WORLDCOVER
    lines = ["# Exploratory: incidence-detrended correlation (NOT pre-registered)\n",
             f"Held-out tiles kept by the registered evaluation ({n_tiles}), class-pure 90 m blocks "
             "as in R1-R3. Local incidence from DEM normals and the zero-Doppler look vector. "
             f"Trend = mean dB in {BIN_DEG:g} deg local-incidence bins (at least {MIN_PER_BIN} "
             "blocks), interpolated; blocks outside the binned range are dropped.\n"]
    for direction, rep in report.items():
        incs = [g[0] for g in rep["tile_geometry"].values()]
        lines += [f"## {direction}\n", f"Ellipsoid incidence over the tiles: {min(incs):.1f}-"
                  f"{max(incs):.1f} deg.\n",
                  "| Class | Blocks | r (raw) | r (incidence removed) | Variance explained by incidence: real / sim | "
                  "Slope, dB per 10 deg (20-60 deg): real / sim |", "|---|---|---|---|---|---|"]
        for c, v in rep["classes"].items():
            lines.append(f"| {name.get(int(c), c)} | {v['blocks_used']} | {v['r_raw']:.3f} | "
                         f"{v['r_detrended']:.3f} | {v['r2_incidence_real']:.2f} / "
                         f"{v['r2_incidence_sim']:.2f} | {v['slope_real_db_per_10deg']:+.2f} / "
                         f"{v['slope_sim_db_per_10deg']:+.2f} |")
        lines += ["", "Median dB per 10 deg of local incidence (real / sim):\n",
                  "| Class | " + " | ".join(f"{a}-{a + 10}" for a in range(0, 90, 10)) + " |",
                  "|---|" + "---|" * 9]
        for c, v in rep["classes"].items():
            cells = {int(lo): f"{r:.1f} / {s:.1f}" for lo, _, r, s in v["curve"]}
            lines.append(f"| {name.get(int(c), c)} | " +
                         " | ".join(cells.get(a, "") for a in range(0, 90, 10)) + " |")
        lines.append("")
    return lines


if __name__ == "__main__":
    sys.exit(main())
