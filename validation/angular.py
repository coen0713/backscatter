"""Round-2 metrics A1-A4 and N (docs/prereg_angular.md, section 4).

Per pass and per class (WorldCover 10, 30, 60, 100), on class-pure 90 m blocks
with no layover or shadow cell in the simulation or the SNAP mask:

  A1  slope of simulated dB against local incidence (20-60 deg), divided by the
      real slope
  A2  median (simulated - real) dB, local incidence below 30 deg
  A3  median (simulated - real) dB, local incidence 40-60 deg
  A4  Pearson r of simulated against real dB, each with its own local-incidence
      trend removed (mean in 2 deg bins of at least 30 blocks, interpolated;
      blocks outside the binned range dropped), as in the exploratory analysis
  N   the round-1 metrics (validation.pipeline.compare / evaluate.pooled)

95% intervals from 1,000 bootstrap resamples of tiles (seed 20261013).
Local incidence: validation.explore_incidence (DEM normal and the zero-Doppler
look vector at the tile centre).

  python -m validation.angular run --region graubuenden \\
      --asc snap_asc.tif --asc-time ... --desc snap_desc.tif --desc-time ... \\
      --table validation/materials_fitted.txt --build build/release --out build/validation/angular

`--region bernese` runs the same metrics on round-1 data (all tiles are
development data there; label the output accordingly).
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "data"))

from fetch_orbit import fetch_orbit  # noqa: E402
from provenance import exe, header  # noqa: E402

from validation import evaluate, explore_incidence, fit, metrics, pipeline, prepare, regions  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
CLASSES = (10, 30, 60, 100)
SLOPE_BAND = (20.0, 60.0)      # A1
LOW_MAX = 30.0                 # A2: local incidence below this
MID_BAND = (40.0, 60.0)        # A3
MIN_BLOCKS = 30                # fewer blocks in a band: metric not computed (NaN)
BOOT_N = 1000
BOOT_SEED = 20261013

# Part A predictions and Part B thresholds (section 4).
PRED_A = {"A1": (1.3, 2.5), "A2_min": 2.0, "A3_abs": 2.0}
THRESH_B = {"A1": (0.75, 1.33), "A2_abs": 2.0, "A3_abs": 1.5}


# ------------------------------------------------------------ blocks

def tile_blocks(sim: pipeline.TileProducts, ref: dict,
                theta: np.ndarray) -> dict[str, np.ndarray]:
    """Class-pure blocks of the four classes, without any layover or shadow cell
    (simulated or SNAP), as dB values plus block-mean local incidence."""
    valid = sim.valid & ref["ref_valid"]
    k = pipeline.BLOCK
    real_lin, ok = metrics.block_reduce(ref["ref_beta0"], valid, k)
    sim_lin, _ = metrics.block_reduce(sim.beta0, valid, k)
    th, _ = metrics.block_reduce(np.nan_to_num(theta, nan=0.0), valid, k)
    theta_ok, _ = metrics.block_reduce(np.isfinite(theta).astype(float), valid, k)
    masked = sim.layover | sim.shadow | ref["ref_layover"] | ref["ref_shadow"]
    clean = metrics.block_fraction(masked, valid, k) == 0.0
    classes = np.zeros(real_lin.shape, dtype=np.int32)
    for c in CLASSES:
        frac = metrics.block_fraction(sim.worldcover == c, valid, k)
        classes[frac * k * k >= pipeline.CLASS_PURITY - 1e-9] = c
    sel = ok & clean & (classes > 0) & (theta_ok == 1.0)
    return {"sim": metrics.to_db(sim_lin)[sel], "real": metrics.to_db(real_lin)[sel],
            "cls": classes[sel], "theta": th[sel]}


# ------------------------------------------------------------ metrics

def _slope(theta: np.ndarray, y: np.ndarray) -> float:
    m = (theta >= SLOPE_BAND[0]) & (theta <= SLOPE_BAND[1])
    if m.sum() < MIN_BLOCKS or np.ptp(theta[m]) == 0:
        return np.nan
    return float(np.polyfit(theta[m], y[m], 1)[0])


def class_metrics(d: dict[str, np.ndarray]) -> dict[str, float]:
    """A1-A4 for one class's blocks (arrays sim, real, theta)."""
    th, s, r = d["theta"], d["sim"], d["real"]
    slope_s, slope_r = _slope(th, s), _slope(th, r)
    low = th < LOW_MAX
    mid = (th >= MID_BAND[0]) & (th <= MID_BAND[1])
    edges = np.arange(0.0, 90.0 + explore_incidence.BIN_DEG, explore_incidence.BIN_DEG)
    ts, tr = explore_incidence.trend(th, s, edges), explore_incidence.trend(th, r, edges)
    ok = np.isfinite(ts) & np.isfinite(tr)
    rs, rr = s[ok] - ts[ok], r[ok] - tr[ok]
    a4 = (float(np.corrcoef(rs, rr)[0, 1]) if ok.sum() >= MIN_BLOCKS and np.std(rs) > 0
          and np.std(rr) > 0 else np.nan)
    return {
        "blocks": int(s.size),
        "slope_sim_db_per_10deg": slope_s * 10, "slope_real_db_per_10deg": slope_r * 10,
        "A1": slope_s / slope_r if np.isfinite(slope_r) and slope_r != 0 else np.nan,
        "A2": float(np.median(s[low] - r[low])) if low.sum() >= MIN_BLOCKS else np.nan,
        "A2_blocks": int(low.sum()),
        "A3": float(np.median(s[mid] - r[mid])) if mid.sum() >= MIN_BLOCKS else np.nan,
        "A3_blocks": int(mid.sum()),
        "A4": a4,
    }


def split_by_class(tiles: list[dict[str, np.ndarray]]) -> dict[int, list[dict[str, np.ndarray]]]:
    """Per class, the per-tile block arrays (for pooling and the tile bootstrap)."""
    out: dict[int, list[dict[str, np.ndarray]]] = {}
    for c in CLASSES:
        out[c] = [{k: v[t["cls"] == c] for k, v in t.items() if k != "cls"} for t in tiles]
    return out


def pool(parts: list[dict[str, np.ndarray]]) -> dict[str, np.ndarray]:
    return {k: np.concatenate([p[k] for p in parts]) for k in ("sim", "real", "theta")}


def bootstrap(parts: list[dict[str, np.ndarray]], n: int = BOOT_N,
              seed: int = BOOT_SEED) -> dict[str, tuple[float, float]]:
    """95% interval of each metric over resamples of tiles (with replacement)."""
    rng = np.random.default_rng(seed)
    draws: dict[str, list[float]] = {k: [] for k in ("A1", "A2", "A3", "A4")}
    for _ in range(n):
        idx = rng.integers(0, len(parts), len(parts))
        m = class_metrics(pool([parts[i] for i in idx]))
        for k in draws:
            draws[k].append(m[k])
    out = {}
    for k, v in draws.items():
        a = np.asarray(v, float)
        a = a[np.isfinite(a)]
        out[k] = ((float(np.percentile(a, 2.5)), float(np.percentile(a, 97.5)))
                  if a.size >= n // 2 else (np.nan, np.nan))
    return out


def pass_metrics(tiles: list[dict[str, np.ndarray]], n_boot: int = BOOT_N) -> dict[int, dict]:
    out = {}
    for c, parts in split_by_class(tiles).items():
        m = class_metrics(pool(parts))
        m["ci95"] = bootstrap(parts, n_boot) if n_boot else {}
        out[c] = m
    return out


# ------------------------------------------------------------ verdicts

def part_a(per_pass: dict[str, dict[int, dict]]) -> dict:
    """Part A: replicates if, in every pass, A1 is within the predicted range for
    at least three of the four classes and A2 >= +2 dB for at least three."""
    detail = {}
    for direction, res in per_pass.items():
        a1 = sum(PRED_A["A1"][0] <= res[c]["A1"] <= PRED_A["A1"][1] for c in CLASSES)
        a2 = sum(res[c]["A2"] >= PRED_A["A2_min"] for c in CLASSES)
        a3 = sum(abs(res[c]["A3"]) <= PRED_A["A3_abs"] for c in CLASSES)
        detail[direction] = {"A1_classes_in_range": int(a1), "A2_classes_at_least_2dB": int(a2),
                             "A3_classes_within_2dB": int(a3),
                             "A1_holds": a1 >= 3, "A2_holds": a2 >= 3, "A3_holds": a3 == len(CLASSES)}
    replicates = all(d["A1_holds"] and d["A2_holds"] for d in detail.values())
    return {"replicates": bool(replicates), "per_pass": detail}


def part_b(per_pass: dict[str, dict[int, dict]], nonreg: dict[str, dict[str, bool]]) -> dict:
    """Part B: succeeds if, in every pass, every class meets the A1-A3 thresholds
    and every round-1 criterion passes."""
    detail = {}
    for direction, res in per_pass.items():
        cls = {c: {"A1": THRESH_B["A1"][0] <= res[c]["A1"] <= THRESH_B["A1"][1],
                   "A2": abs(res[c]["A2"]) <= THRESH_B["A2_abs"],
                   "A3": abs(res[c]["A3"]) <= THRESH_B["A3_abs"]} for c in CLASSES}
        detail[direction] = {"classes": cls, "nonregression": nonreg[direction]}
    ok = all(all(all(v.values()) for v in d["classes"].values()) and all(d["nonregression"].values())
             for d in detail.values())
    return {"succeeds": bool(ok), "per_pass": detail}


# ------------------------------------------------------------ driver

def run(args: argparse.Namespace) -> int:
    tiles = regions.region_tiles(args.region)
    if args.tiles:
        tiles = [t for t in tiles if t.name in args.tiles]
    for t in tiles:
        pipeline.check_allowed(t)
    render = exe(args.build, "apps", "bsar_render")
    fit.read_table(args.table)  # fails early on a malformed table
    args.out.mkdir(parents=True, exist_ok=True)
    blocks: dict[str, dict[str, dict]] = {}
    compared: dict[str, list[dict]] = {}
    for direction, snap, utc in (("ascending", args.asc, args.asc_time),
                                 ("descending", args.desc, args.desc_time)):
        when = datetime.fromisoformat(utc)
        orbit = fetch_orbit(args.mission, when, ROOT / "data" / "raw")
        ts, sv, t0 = explore_incidence.read_eof(orbit)
        for tile in tiles:
            prep = pipeline.prepare_tile(tile, args.out / direction / tile.name)
            sim = pipeline.render_tile(prep, orbit, utc, render, material_table=args.table)
            ref = evaluate.load_snap_reference(snap, prep.grid, args.beta0_band, args.mask_band,
                                               args.layover_bit, args.shadow_bit)
            w, s, e, n = tile.bounds
            lat, lon = 0.5 * (s + n), 0.5 * (w + e)
            h = float(np.nanmedian(np.where(prep.dem == prepare.NODATA, np.nan, prep.dem)))
            inc, az = explore_incidence.look_geometry(ts, sv, (when - t0).total_seconds(), lat, lon, h)
            theta = explore_incidence.local_incidence(prep.dem, prep.grid, inc, az, lat, lon)
            blocks.setdefault(direction, {})[tile.name] = tile_blocks(sim, ref, theta)
            compared.setdefault(direction, []).append(pipeline.compare(sim, **ref))
            print(f"{direction} {tile.name}: incidence {inc:.2f}", flush=True)
    # Tile exclusions as round 1 with amendment 1: failing in either pass excludes from both.
    combined = pipeline.combine_exclusions(
        {d: {r["tile"]: r["exclusions"] for r in rs} for d, rs in compared.items()})
    excluded = {t: r for t, r in combined.items() if r}
    per_pass, nonreg, pooled_n = {}, {}, {}
    for direction in blocks:
        kept = [b for name, b in blocks[direction].items() if name not in excluded]
        per_pass[direction] = pass_metrics(kept, args.bootstrap)
        for r in compared[direction]:
            r["exclusions"] = combined[r["tile"]]
        p = evaluate.pooled(compared[direction])
        pooled_n[direction] = p
        nonreg[direction] = {k: bool(v) for k, v in evaluate.verdicts(p).items()}
    result = {"region": args.region, "table": str(args.table), "excluded": excluded,
              "metrics": per_pass, "part_a": part_a(per_pass), "part_b": part_b(per_pass, nonreg),
              "nonregression_pooled": pooled_n}
    (args.out / "results.json").write_text(json.dumps(result, indent=1, default=str))
    lines = report_lines(result, header("python -m validation.angular", args.build), args.label)
    (args.out / "report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    return 0


def _ci(m: dict, k: str, fmt: str) -> str:
    v = m[k]
    lo, hi = m.get("ci95", {}).get(k, (np.nan, np.nan))
    s = format(v, fmt) if np.isfinite(v) else "n/a"
    return s + (f" [{format(lo, fmt)}, {format(hi, fmt)}]" if np.isfinite(lo) else "")


def report_lines(result: dict, head: str, label: str) -> list[str]:
    names = metrics.WORLDCOVER
    lines = [f"# Angular-response metrics: {result['region']} ({label})\n", head,
             f"Material table: `{result['table']}`.\n"]
    if result["excluded"]:
        lines.append(f"Excluded tiles (section 3, either pass): {len(result['excluded'])}: "
                     + ", ".join(sorted(result["excluded"])) + "\n")
    for direction, res in result["metrics"].items():
        lines += [f"## {direction}\n",
                  "| Class | Blocks | A1 slope ratio | Slope real / sim (dB per 10°) | A2 bias < 30° (blocks) "
                  "| A3 bias 40-60° (blocks) | A4 r, incidence removed |", "|---|---|---|---|---|---|---|"]
        for c in CLASSES:
            m = res[c]
            lines.append(f"| {names[c]} | {m['blocks']} | {_ci(m, 'A1', '.2f')} | "
                         f"{m['slope_real_db_per_10deg']:+.2f} / {m['slope_sim_db_per_10deg']:+.2f} | "
                         f"{_ci(m, 'A2', '+.2f')} dB ({m['A2_blocks']}) | "
                         f"{_ci(m, 'A3', '+.2f')} dB ({m['A3_blocks']}) | {_ci(m, 'A4', '.3f')} |")
        n = result["part_b"]["per_pass"][direction]["nonregression"]
        lines += ["", "Round-1 criteria on these tiles: " +
                  ", ".join(f"{k} {'pass' if v else 'FAIL'}" for k, v in n.items()), ""]
    a = result["part_a"]
    lines += ["## Verdicts\n",
              f"Part A (replication; A1 in 1.3-2.5 and A2 >= +2 dB for at least 3 of 4 classes "
              f"in both passes): **{'replicates' if a['replicates'] else 'does not replicate'}**.",
              "Part A detail: " + "; ".join(
                  f"{d}: A1 {v['A1_classes_in_range']}/4, A2 {v['A2_classes_at_least_2dB']}/4, "
                  f"A3 within 2 dB {v['A3_classes_within_2dB']}/4" for d, v in a["per_pass"].items()),
              f"Part B thresholds (only meaningful for a revised model): "
              f"**{'met' if result['part_b']['succeeds'] else 'not met'}**.", ""]
    return lines


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="command", required=True)
    r = sub.add_parser("run")
    r.add_argument("--region", choices=sorted(regions.REGIONS), required=True)
    r.add_argument("--tiles", nargs="*", help="restrict to these tile names")
    r.add_argument("--asc", type=Path, required=True)
    r.add_argument("--asc-time", required=True)
    r.add_argument("--desc", type=Path, required=True)
    r.add_argument("--desc-time", required=True)
    r.add_argument("--mission", default="S1A")
    r.add_argument("--beta0-band", default="1")
    r.add_argument("--mask-band", default="2")
    r.add_argument("--layover-bit", type=int, default=1)
    r.add_argument("--shadow-bit", type=int, default=2)
    r.add_argument("--table", type=Path, default=ROOT / "validation" / "materials_fitted.txt")
    r.add_argument("--build", type=Path, default=ROOT / "build" / "release")
    r.add_argument("--bootstrap", type=int, default=BOOT_N, help="resamples (0 = skip)")
    r.add_argument("--label", default="test", help="e.g. 'development, round-1 data'")
    r.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    return run(args)


if __name__ == "__main__":
    sys.exit(main())
