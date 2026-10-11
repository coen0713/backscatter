"""Per-class fitting of the diffuse parameters k and n (docs/prereg.md, section 5).

Only `k` (diffuse coefficient) and `n` (diffuse exponent) of each WorldCover
class may be fitted, and only on calibration tiles. Everything else in the
material table stays fixed.

Method. The geometric renderer is linear in the diffuse coefficients: ray paths,
specular reflections and bounce weights depend only on permittivity and RMS
slope. A block's signal is split into

  S1    = single-bounce render with k = 0 for every class (specular lobe)
  D1_n  = single-bounce render with k = 1, exponent n, minus S1, for n on N_GRID
  M     = (2-bounce render) - (1-bounce render) with the *current* table:
          everything that bounced more than once

and a 90 m block of class c is predicted as  S1 + k_c * D1_{n_c} + M.  Single
bounce is attributed to the block's own class; multi-bounce is not (a lake's
multi-bounce signal is mostly sensor -> water -> cliff -> sensor and scales
with the cliff's k), so M is held fixed during a fit and re-rendered with the
fitted table for the next iteration (ITERATIONS rounds). Between grid values of
n, log D1 is interpolated linearly in n (exact for a block with a single local
incidence angle, since D1 ~ cos^n).

Blocks used: interior blocks (the class on all 9 cells and on a 1-cell ring
around them), every cell valid in simulation and
reference, no layover and no shadow cell (layover mixes classes; shadow has no
signal). Snow and ice (70) are excluded, as in the radiometry metrics.

Objective (pre-registered): the median absolute dB error over those blocks,
pooled over all calibration tiles and both passes; minimised over k for each
n, then over n.

  python -m validation.fit --asc snap_asc.tif --asc-time 2024-08-07T17:24:01 \\
      --desc snap_desc.tif --desc-time 2024-08-06T05:35:49 --build build/release \\
      --out build/validation/fit
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass, replace
from datetime import datetime
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "data"))

from fetch_orbit import fetch_orbit  # noqa: E402
from provenance import exe, header  # noqa: E402

from validation import metrics, pipeline, regions  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent

N_GRID = (0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0)
MIN_BLOCKS = 100   # fewer pure blocks than this: keep the starting values
SNOW_ICE = 70
INTERIOR_RING = 1  # cells of the same class required around a fitting block
ITERATIONS = 2     # multi-bounce re-render / refit rounds
N_WEAK_DB = 0.1    # objective rise at n +- 0.5 below which n counts as weakly constrained


# ------------------------------------------------------------------ tables

@dataclass(frozen=True)
class MaterialEntry:
    code: int
    name: str
    eps_re: float
    eps_im: float
    k: float
    n: float
    rms_slope: float
    pec: bool = False


def read_table(path: Path) -> list[MaterialEntry]:
    """Read a material table whose entries are all `custom`."""
    entries = []
    for raw in path.read_text().splitlines():
        line = raw.split("#", 1)[0].split()
        if not line:
            continue
        if line[1] != "custom":
            raise ValueError(f"{path}: entry {line[0]} must be 'custom' to be fitted")
        entries.append(MaterialEntry(int(line[0]), line[2], float(line[3]), float(line[4]),
                                     float(line[5]), float(line[6]), float(line[7]),
                                     len(line) > 8 and line[8] == "pec"))
    return entries


def write_table(path: Path, entries: list[MaterialEntry], comment: str) -> None:
    lines = [f"# {c}" for c in comment.splitlines()]
    lines.append("# code  kind    name            eps_re  eps_im  k        n      rms_slope")
    for e in entries:
        lines.append(f"{e.code:<7d} custom  {e.name:<15s} {e.eps_re:<7g} {e.eps_im:<7g} "
                     f"{e.k:<8.5g} {e.n:<6.4g} {e.rms_slope:g}" + (" pec" if e.pec else ""))
    path.write_text("\n".join(lines) + "\n")


def uniform(entries: list[MaterialEntry], k: float, n: float) -> list[MaterialEntry]:
    return [replace(e, k=k, n=n) for e in entries]


# ------------------------------------------------------------------ basis

@dataclass
class Basis:
    """Fitting basis of one tile and pass, on the tile's UTM grid."""

    tile: regions.Tile
    specular: np.ndarray            # S1
    diffuse: dict[float, np.ndarray]  # D1_n
    usable: np.ndarray              # valid, no layover, no shadow
    worldcover: np.ndarray


def _check_calibration(prep: pipeline.PreparedTile) -> None:
    if prep.tile not in regions.calibration_tiles():
        raise PermissionError(f"{prep.tile.name}: fitting uses calibration tiles only")


def render_basis(prep: pipeline.PreparedTile, orbit: Path, utc: str, render: Path,
                 base: list[MaterialEntry], n_grid=N_GRID) -> Basis:
    """S1 and D1_n (single bounce)."""
    _check_calibration(prep)
    tables = prep.work / "fit_tables"
    tables.mkdir(exist_ok=True)
    zero = tables / "k0.txt"
    write_table(zero, uniform(base, 0.0, 1.0), "fit basis: diffuse off")
    spec = pipeline.render_tile(prep, orbit, utc, render, name="fit_specular",
                                material_table=zero, bounces=1)
    diffuse = {}
    for n in n_grid:
        table = tables / f"k1_n{n:g}.txt"
        write_table(table, uniform(base, 1.0, n), f"fit basis: k = 1, n = {n:g}")
        r = pipeline.render_tile(prep, orbit, utc, render, name=f"fit_n{n:g}",
                                 material_table=table, bounces=1)
        diffuse[n] = np.maximum(r.beta0 - spec.beta0, 0.0)
    usable = spec.valid & ~spec.layover & ~spec.shadow
    return Basis(prep.tile, spec.beta0, diffuse, usable, prep.worldcover)


def multi_bounce(prep: pipeline.PreparedTile, orbit: Path, utc: str, render: Path,
                 table: list[MaterialEntry], tag: str) -> np.ndarray:
    """M: the 2-bounce render minus the 1-bounce render with `table`."""
    _check_calibration(prep)
    path = prep.work / "fit_tables" / f"current_{tag}.txt"
    path.parent.mkdir(exist_ok=True)
    write_table(path, table, f"fit: current parameters ({tag})")
    two = pipeline.render_tile(prep, orbit, utc, render, name=f"fit_{tag}_b2",
                               material_table=path, bounces=2)
    one = pipeline.render_tile(prep, orbit, utc, render, name=f"fit_{tag}_b1",
                               material_table=path, bounces=1)
    return two.beta0 - one.beta0


# ------------------------------------------------------------------ blocks

@dataclass
class ClassBlocks:
    """Block means (linear) for one class, pooled over tiles and passes."""

    y: list[np.ndarray]
    s: list[np.ndarray]
    f: dict[float, list[np.ndarray]]

    def arrays(self) -> tuple[np.ndarray, np.ndarray, dict[float, np.ndarray]]:
        return (np.concatenate(self.y), np.concatenate(self.s),
                {n: np.concatenate(v) for n, v in self.f.items()})


def collect_blocks(basis: Basis, ref_beta0: np.ndarray, ref_valid: np.ndarray,
                   pooled: dict[int, ClassBlocks], offset: np.ndarray | None = None) -> None:
    """Add this tile's fitting blocks to `pooled`. `offset` (M, the multi-bounce
    part) is added to the fixed part of the prediction."""
    usable = basis.usable & ref_valid
    k = pipeline.BLOCK
    fixed = basis.specular if offset is None else basis.specular + offset
    y_blk, ok = metrics.block_reduce(ref_beta0, usable, k, min_fraction=1.0)
    s_blk, _ = metrics.block_reduce(fixed, usable, k, min_fraction=1.0)
    f_blk = {n: metrics.block_reduce(f, usable, k, min_fraction=1.0)[0]
             for n, f in basis.diffuse.items()}
    for c in np.unique(basis.worldcover):
        if c in (0, SNOW_ICE):
            continue
        # Interior blocks only: the class on all 9 cells and on a 1-cell ring
        # around them. With 7-of-9 purity, two bright neighbour cells (or one
        # cell of resampling bleed at a class edge) biased dark classes upwards.
        interior = ~metrics.dilate(basis.worldcover != c, INTERIOR_RING)
        pure = metrics.block_fraction(interior, np.ones_like(usable), k)
        sel = ok & (pure >= 1.0 - 1e-9) & (y_blk > 0)
        if not sel.any():
            continue
        entry = pooled.setdefault(int(c), ClassBlocks([], [], {n: [] for n in basis.diffuse}))
        entry.y.append(y_blk[sel])
        entry.s.append(s_blk[sel])
        for n in basis.diffuse:
            entry.f[n].append(f_blk[n][sel])


# ------------------------------------------------------------------ fitting

def diffuse_at(f: dict[float, np.ndarray], n: float) -> np.ndarray:
    """F_n for any n in the grid's range: log-linear interpolation between nodes."""
    grid = sorted(f)
    if n <= grid[0]:
        return f[grid[0]]
    if n >= grid[-1]:
        return f[grid[-1]]
    hi = next(g for g in grid if g >= n)
    lo = grid[grid.index(hi) - 1]
    if hi == n:
        return f[hi]
    w = (n - lo) / (hi - lo)
    a = np.log(np.maximum(f[lo], 1e-30))
    b = np.log(np.maximum(f[hi], 1e-30))
    return np.exp((1 - w) * a + w * b)


def objective(y: np.ndarray, s: np.ndarray, fn: np.ndarray, k: float) -> float:
    """Median absolute dB error of the prediction s + k fn against y."""
    return float(np.median(np.abs(metrics.to_db(s + k * fn) - metrics.to_db(y))))


def _golden(fun, lo: float, hi: float, tol: float) -> float:
    g = (np.sqrt(5.0) - 1.0) / 2.0
    a, b = lo, hi
    c, d = b - g * (b - a), a + g * (b - a)
    fc, fd = fun(c), fun(d)
    while b - a > tol:
        if fc < fd:
            b, d, fd = d, c, fc
            c = b - g * (b - a)
            fc = fun(c)
        else:
            a, c, fc = c, d, fd
            d = a + g * (b - a)
            fd = fun(d)
    return 0.5 * (a + b)


def best_k(y: np.ndarray, s: np.ndarray, fn: np.ndarray) -> tuple[float, float]:
    """(k, objective): coarse log-spaced scan, then golden-section in log k."""
    ks = np.geomspace(1e-5, 10.0, 121)
    vals = [objective(y, s, fn, k) for k in ks]
    i = int(np.argmin(vals))
    lo, hi = np.log(ks[max(i - 1, 0)]), np.log(ks[min(i + 1, len(ks) - 1)])
    logk = _golden(lambda lk: objective(y, s, fn, float(np.exp(lk))), lo, hi, 1e-4)
    k = float(np.exp(logk))
    return k, objective(y, s, fn, k)


def fit_class(y: np.ndarray, s: np.ndarray, f: dict[float, np.ndarray]) -> dict:
    """Fit (k, n) for one class: k by 1-D search for each n, n by a scan at
    0.05 steps over the grid's range then golden-section refinement."""
    grid = sorted(f)
    ns = np.arange(grid[0], grid[-1] + 1e-9, 0.05)
    scan = [best_k(y, s, diffuse_at(f, float(n)))[1] for n in ns]
    i = int(np.argmin(scan))
    lo, hi = ns[max(i - 1, 0)], ns[min(i + 1, len(ns) - 1)]
    n = _golden(lambda v: best_k(y, s, diffuse_at(f, v))[1], float(lo), float(hi), 1e-3)
    k, obj = best_k(y, s, diffuse_at(f, n))
    # How well is n constrained? A class seen over a narrow range of local
    # incidence (a flat lake) cannot separate k from n: k cos^n is nearly the
    # same for many (k, n) pairs. Report the objective's rise at n +- 0.5
    # (with k re-optimised); below N_WEAK_DB, n is weakly constrained.
    rises = [best_k(y, s, diffuse_at(f, float(np.clip(n + d, grid[0], grid[-1]))))[1] - obj
             for d in (-0.5, 0.5) if grid[0] <= n + d <= grid[-1]]
    rise = float(min(rises)) if rises else float("nan")
    return {"k": k, "n": float(n), "median_abs_db": obj, "blocks": int(y.size),
            "n_profile_rise_db": rise, "n_weak": bool(rise < N_WEAK_DB)}


def fit_all(pooled: dict[int, ClassBlocks], base: list[MaterialEntry]) -> tuple[list[MaterialEntry], dict]:
    by_code = {e.code: e for e in base}
    fitted = list(base)
    report = {}
    for code, blocks in sorted(pooled.items()):
        y, s, f = blocks.arrays()
        start = by_code.get(code)
        if start is None:
            continue
        before = objective(y, s, diffuse_at(f, start.n), start.k)
        if y.size < MIN_BLOCKS:
            report[code] = {"name": start.name, "blocks": int(y.size), "kept_start": True,
                            "k": start.k, "n": start.n, "median_abs_db": before,
                            "start_median_abs_db": before}
            continue
        res = fit_class(y, s, f)
        res.update(name=start.name, kept_start=False, start_k=start.k, start_n=start.n,
                   start_median_abs_db=before)
        report[code] = res
        fitted = [replace(e, k=res["k"], n=res["n"]) if e.code == code else e for e in fitted]
    return fitted, report


@dataclass
class FitContext:
    """One tile and pass: its basis and the reference to fit against."""

    prep: pipeline.PreparedTile
    orbit: Path
    utc: str
    basis: Basis
    ref_beta0: np.ndarray
    ref_valid: np.ndarray


def run_fit(contexts: list[FitContext], base: list[MaterialEntry], render: Path,
            iterations: int = ITERATIONS
            ) -> tuple[list[MaterialEntry], dict, list[dict], dict[int, ClassBlocks]]:
    """Fit k and n per class. Each iteration re-renders the multi-bounce part M
    with the current table, then refits every class with M held fixed.
    Returns (fitted table, final report, per-iteration reports, final blocks)."""
    table = list(base)
    history = []
    pooled: dict[int, ClassBlocks] = {}
    report: dict = {}
    for it in range(iterations):
        pooled = {}
        for ctx in contexts:
            offset = multi_bounce(ctx.prep, ctx.orbit, ctx.utc, render, table, f"it{it}")
            collect_blocks(ctx.basis, ctx.ref_beta0, ctx.ref_valid, pooled, offset)
        # Starting values stay the reference for "kept" classes and for the
        # before/after comparison; the fit itself starts from scratch each round.
        table, report = fit_all(pooled, base)
        history.append(report)
    return table, report, history, pooled


def report_lines(report: dict) -> list[str]:
    lines = ["| Class | Blocks | Start k / n | Fitted k / n | Median abs. error, start -> fitted "
             "| n constrained? (objective rise at n +- 0.5) |",
             "|---|---|---|---|---|---|"]
    for code, r in sorted(report.items()):
        if r.get("kept_start"):
            lines.append(f"| {code} {r['name']} | {r['blocks']} | {r['k']:.4g} / {r['n']:.2f} | "
                         f"kept (fewer than {MIN_BLOCKS} blocks) | {r['median_abs_db']:.2f} dB | |")
        else:
            lines.append(f"| {code} {r['name']} | {r['blocks']} | {r['start_k']:.4g} / "
                         f"{r['start_n']:.2f} | {r['k']:.4g} / {r['n']:.2f} | "
                         f"{r['start_median_abs_db']:.2f} -> {r['median_abs_db']:.2f} dB | "
                         f"{'weakly' if r['n_weak'] else 'yes'} ({r['n_profile_rise_db']:.2f} dB) |")
    lines += ["", "Where n is only weakly constrained (the class is seen over a narrow range of "
              "local incidence), k and n trade off and only their combination k cos^n is "
              "determined; the fitted pair is still the best predictor over that range, but n "
              "should not be read as a physical result."]
    return lines


# ------------------------------------------------------------------ CLI

def main() -> int:
    from validation import evaluate  # noqa: PLC0415 (heavy import only for the CLI)

    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--tiles", nargs="*", help="restrict to these calibration tile names")
    p.add_argument("--asc", type=Path)
    p.add_argument("--asc-time")
    p.add_argument("--desc", type=Path)
    p.add_argument("--desc-time")
    p.add_argument("--mission", default="S1A")
    p.add_argument("--beta0-band", default="Beta0_VV")
    p.add_argument("--mask-band", default="layover_shadow_mask")
    p.add_argument("--layover-bit", type=int, default=1)
    p.add_argument("--shadow-bit", type=int, default=2)
    p.add_argument("--table", type=Path, default=pipeline.MATERIAL_TABLE)
    p.add_argument("--build", type=Path, default=ROOT / "build" / "release")
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()

    tiles = regions.calibration_tiles()
    if args.tiles:
        unknown = set(args.tiles) - {t.name for t in tiles}
        if unknown:
            raise SystemExit(f"not calibration tiles (fitting uses calibration only): {sorted(unknown)}")
        tiles = [t for t in tiles if t.name in args.tiles]
    render = exe(args.build, "apps", "bsar_render")
    base = read_table(args.table)
    contexts = []
    for direction, snap, utc in (("ascending", args.asc, args.asc_time),
                                 ("descending", args.desc, args.desc_time)):
        if snap is None:
            continue
        orbit = fetch_orbit(args.mission, datetime.fromisoformat(utc), ROOT / "data" / "raw")
        for tile in tiles:
            prep = pipeline.prepare_tile(tile, args.out / direction / tile.name)
            basis = render_basis(prep, orbit, utc, render, base)
            ref = evaluate.load_snap_reference(snap, prep.grid, args.beta0_band, args.mask_band,
                                               args.layover_bit, args.shadow_bit)
            contexts.append(FitContext(prep, orbit, utc, basis, ref["ref_beta0"], ref["ref_valid"]))
            print(f"basis {direction} {tile.name}", flush=True)
    fitted, report, _, _ = run_fit(contexts, base, render)
    args.out.mkdir(parents=True, exist_ok=True)
    write_table(args.out / "materials_fitted.txt", fitted,
                "Fitted on calibration tiles by python -m validation.fit (k and n only).")
    (args.out / "fit.json").write_text(json.dumps(report, indent=1))
    lines = ["# Per-class fit (calibration tiles)\n", header("python -m validation.fit", args.build)]
    lines += report_lines(report)
    (args.out / "report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    sys.exit(main())
