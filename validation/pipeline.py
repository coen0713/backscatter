"""Per-tile simulation and comparison (docs/prereg.md, sections 4 and 6).

simulate_tile: DEM crop with margin -> bsar_render (real orbit, --geocode) ->
products resampled to the tile's UTM 30 m grid, plus WorldCover and the DEM.

compare: the pre-registered metrics between a simulated tile and a reference
(real SNAP output, or a synthetic stand-in during the dry run).

Held-out tiles are refused unless docs/prereg.md is marked REGISTERED, and
round-2 test tiles unless docs/prereg_angular.md is.
"""

from __future__ import annotations

import re
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from rasterio.enums import Resampling

from validation import metrics, prepare, regions

ROOT = Path(__file__).resolve().parent.parent
PREREG = ROOT / "docs" / "prereg.md"
PREREGS = {"bernese": PREREG, "graubuenden": ROOT / "docs" / "prereg_angular.md",
           "tyrol": ROOT / "docs" / "prereg_partb.md"}
MATERIAL_TABLE = Path(__file__).resolve().parent / "materials_worldcover.txt"

MARGIN_DEG = 0.05      # terrain outside the tile that can shadow or lay over into it
RELIEF_MIN_M = 1000.0  # "mountainous" tile threshold for the geometry metrics
BLOCK = 3              # 3 x 3 cells of 30 m = 90 m blocks
CLASS_PURITY = 7       # a block belongs to a class if >= 7 of its 9 cells do
SNOW_ICE = 70
WATER = 80


def prereg_registered(path: Path = PREREG) -> bool:
    text = path.read_text(encoding="utf-8")
    status = re.search(r"\*\*Status:\s*([^*]+)\*\*", text)
    return bool(status) and status.group(1).strip().upper().startswith("REGISTERED")


def check_allowed(tile: regions.Tile) -> None:
    if not tile.in_region:
        raise ValueError(f"tile ({tile.i}, {tile.j}) lies outside region {tile.region}")
    if tile.region != "bernese":
        # Round 2: every tile is a test tile, refused until its prereg is registered.
        if not prereg_registered(PREREGS[tile.region]):
            raise PermissionError(
                f"{tile.name} is a test tile of region {tile.region}; "
                f"{PREREGS[tile.region].name} is not registered yet.")
        return
    if tile in regions.calibration_tiles():
        return
    if not prereg_registered():
        raise PermissionError(
            f"{tile.name} is a held-out tile; docs/prereg.md is not registered yet. "
            "Only calibration tiles may be processed before registration.")


@dataclass
class TileProducts:
    """Everything on one tile's UTM 30 m grid. Masks are bool arrays."""

    tile: regions.Tile
    grid: prepare.Grid
    beta0: np.ndarray
    layover: np.ndarray
    shadow: np.ndarray
    valid: np.ndarray
    dem: np.ndarray
    worldcover: np.ndarray
    info: dict = field(default_factory=dict)


@dataclass
class PreparedTile:
    """Inputs for rendering one tile, independent of orbit and materials."""

    tile: regions.Tile
    work: Path
    dem_asc: Path
    classes_asc: Path
    grid: prepare.Grid  # the tile's UTM 30 m grid
    dem: np.ndarray     # ellipsoidal heights on `grid`
    worldcover: np.ndarray  # class codes on `grid`


def prepare_tile(tile: regions.Tile, work: Path) -> PreparedTile:
    """DEM crop with margin, class grid on the DEM grid, and the tile's UTM grid."""
    check_allowed(tile)
    work.mkdir(parents=True, exist_ok=True)
    w, s, e, n = tile.bounds
    heights, dem_grid = prepare.ellipsoidal_dem((w - MARGIN_DEG, s - MARGIN_DEG,
                                                 e + MARGIN_DEG, n + MARGIN_DEG))
    dem_asc = work / "dem.asc"
    prepare.write_asc(dem_asc, heights, dem_grid)
    classes_asc = work / "classes.asc"
    prepare.write_asc(classes_asc, prepare.worldcover(dem_grid).astype(np.float64), dem_grid)
    grid = prepare.utm_grid(tile.bounds)
    dem_utm = prepare.resample(heights, dem_grid, grid, Resampling.bilinear, prepare.NODATA)
    return PreparedTile(tile, work, dem_asc, classes_asc, grid, dem_utm, prepare.worldcover(grid))


def render_tile(prep: PreparedTile, orbit: Path, utc: str, render: Path, name: str = "sim",
                material_table: Path = MATERIAL_TABLE, range_spacing: float = 10.0,
                azimuth_spacing: float = 10.0, bounces: int = 2,
                threads: int = 0) -> TileProducts:
    """bsar_render with the real orbit and a material table, geocoded and
    resampled to the tile's UTM grid. `name` keeps several renders of one
    tile (e.g. the fitting basis) apart."""
    check_allowed(prep.tile)
    prefix = prep.work / name
    cmd = [str(render), "--scene", f"dem:{prep.dem_asc}", "--geographic", "--classes",
           str(prep.classes_asc), "--material-table", str(material_table), "--orbit", str(orbit),
           "--time", utc, "--range-spacing", str(range_spacing), "--azimuth-spacing",
           str(azimuth_spacing), "--bounces", str(bounces), "--threads", str(threads),
           "--geocode", "--out", str(prefix)]
    # Resumable: skip the render if this exact command already completed (the
    # stamp records the command and the table's contents).
    stamp = prep.work / f"{name}.done"
    signature = " ".join(cmd) + "\n" + Path(material_table).read_text()
    log = prep.work / f"{name}.log"
    if not (stamp.exists() and stamp.read_text() == signature and log.exists()):
        run = subprocess.run(cmd, capture_output=True, text=True)
        log.write_text(run.stdout + run.stderr)
        if run.returncode != 0:
            raise RuntimeError(f"bsar_render failed for {prep.tile.name}:\n{run.stderr}")
        stamp.write_text(signature)

    out = {}
    for product, resampling in (("beta0", Resampling.bilinear), ("layover", Resampling.nearest),
                                ("shadow", Resampling.nearest), ("valid", Resampling.nearest)):
        values, src_grid, nodata = prepare.read_raster(Path(f"{prefix}_geo_{product}.asc"))
        out[product] = prepare.resample(values, src_grid, prep.grid, resampling, nodata)
    valid = ((out["valid"] > 0.5) & (out["beta0"] != prepare.NODATA) &
             (prep.dem != prepare.NODATA))
    return TileProducts(
        tile=prep.tile, grid=prep.grid, beta0=np.where(valid, out["beta0"], 0.0),
        layover=(out["layover"] > 0.5) & valid, shadow=(out["shadow"] > 0.5) & valid,
        valid=valid, dem=prep.dem, worldcover=prep.worldcover,
        info={"render_log": log.read_text().strip().splitlines()})


def simulate_tile(tile: regions.Tile, orbit: Path, utc: str, render: Path, work: Path,
                  material_table: Path = MATERIAL_TABLE, **kwargs) -> TileProducts:
    return render_tile(prepare_tile(tile, work), orbit, utc, render,
                       material_table=material_table, **kwargs)


def exclusion_reasons(dem: np.ndarray, worldcover: np.ndarray, sim_valid: np.ndarray,
                      ref_valid: np.ndarray) -> list[str]:
    """Pre-registered tile-level exclusions (section 3) for one pass. Empty list = keep."""
    reasons = []
    dem_nodata = np.mean(dem == prepare.NODATA)
    if dem_nodata > 0.05:
        reasons.append(f"DEM no-data {dem_nodata:.1%} > 5%")
    water = np.mean(worldcover == WATER)
    if water > 0.5:
        reasons.append(f"water {water:.1%} > 50%")
    if not (sim_valid & ref_valid).all():
        reasons.append("not fully inside both footprints")
    return reasons


def tile_exclusions(sim: TileProducts, ref_valid: np.ndarray) -> list[str]:
    return exclusion_reasons(sim.dem, sim.worldcover, sim.valid, ref_valid)


def combine_exclusions(per_pass: dict[str, dict[str, list[str]]]) -> dict[str, list[str]]:
    """Section 3 excludes a tile "entirely" if it fails in either scene, so the
    per-pass reasons {direction: {tile: reasons}} become {tile: reasons} for
    every pass, each reason tagged with the pass it came from."""
    out: dict[str, list[str]] = {}
    for direction, tiles in per_pass.items():
        for tile, reasons in tiles.items():
            out.setdefault(tile, []).extend(f"{direction}: {r}" for r in reasons)
    return out


def compare(sim: TileProducts, ref_beta0: np.ndarray, ref_layover: np.ndarray,
            ref_shadow: np.ndarray, ref_valid: np.ndarray) -> dict:
    """The pre-registered metrics for one tile (section 4). Linear beta0 in,
    dB inside. Returns a flat dict; NaN where a metric does not apply."""
    valid = sim.valid & ref_valid
    dem_ok = sim.dem != prepare.NODATA
    relief = float(sim.dem[dem_ok].max() - sim.dem[dem_ok].min()) if dem_ok.any() else 0.0
    mountainous = relief >= RELIEF_MIN_M
    res: dict = {"tile": sim.tile.name, "relief_m": relief, "mountainous": mountainous,
                 "valid_fraction": float(valid.mean()),
                 "exclusions": tile_exclusions(sim, ref_valid)}

    # G1, G2: masks against the reference, mountainous tiles only.
    res["G1_layover_iou"] = metrics.iou(sim.layover, ref_layover, valid) if mountainous else np.nan
    res["G2_shadow_iou"] = metrics.iou(sim.shadow, ref_shadow, valid) if mountainous else np.nan
    # Same masks, tolerating up to one 30 m cell of misregistration.
    res["G1_layover_f1"] = (metrics.tolerant_f1(sim.layover, ref_layover, valid)
                            if mountainous else np.nan)
    res["G2_shadow_f1"] = (metrics.tolerant_f1(sim.shadow, ref_shadow, valid)
                           if mountainous else np.nan)

    # G3: real beta0 in predicted-shadow 90 m blocks (>= 50% shadow) vs blocks
    # with no predicted shadow.
    real_db = metrics.to_db(ref_beta0)
    sim_db = metrics.to_db(sim.beta0)
    # 90 m blocks: average linear power (multilooking), then convert to dB.
    real_lin, ok = metrics.block_reduce(ref_beta0, valid, BLOCK)
    sim_lin, _ = metrics.block_reduce(sim.beta0, valid, BLOCK)
    real_blk = metrics.to_db(real_lin)
    sim_blk = metrics.to_db(sim_lin)
    shadow_frac = metrics.block_fraction(sim.shadow, valid, BLOCK)
    pos = real_blk[ok & (shadow_frac >= 0.5)]
    neg = real_blk[ok & (shadow_frac == 0.0)]
    res["G3_shadow_auc"] = metrics.auc_lower_is_positive(pos, neg) if mountainous else np.nan
    res["G3_blocks"] = (int(pos.size), int(neg.size))

    # R1-R3: per-class correlation on class-pure blocks, snow/ice excluded.
    classes = np.zeros(real_blk.shape, dtype=np.int32)
    for c in np.unique(sim.worldcover):
        if c in (0, SNOW_ICE):
            continue
        frac = metrics.block_fraction(sim.worldcover == c, valid, BLOCK)
        classes[frac * BLOCK * BLOCK >= CLASS_PURITY - 1e-9] = c
    res["R_per_class"] = metrics.per_class_pearson(sim_blk, real_blk, classes, ok & (classes > 0),
                                                   min_count=30)

    # L1: geolocation offset of the real image relative to the simulation.
    fill_s = np.where(valid, sim_db, np.median(sim_db[valid]) if valid.any() else 0.0)
    fill_r = np.where(valid, real_db, np.median(real_db[valid]) if valid.any() else 0.0)
    res["L1_offset_px"] = metrics.phase_correlation_offset(fill_s, fill_r)

    # Raw pieces for pooling across tiles (validation/evaluate.py); not reported.
    res["_arrays"] = {
        "valid": valid, "mountainous": mountainous,
        "sim_layover": sim.layover, "ref_layover": ref_layover.astype(bool),
        "sim_shadow": sim.shadow, "ref_shadow": ref_shadow.astype(bool),
        "g3_pos": pos, "g3_neg": neg,
        "blk_sim": sim_blk[ok & (classes > 0)], "blk_real": real_blk[ok & (classes > 0)],
        "blk_class": classes[ok & (classes > 0)],
    }
    return res
