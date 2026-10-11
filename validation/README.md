# Real-data validation pipeline

Implements the protocol in [`docs/prereg.md`](../docs/prereg.md): simulate each
0.1° tile with a real DEM and a real Sentinel-1 orbit, put simulation and real
image on one UTM 30 m grid, and compute the pre-registered metrics.

**Held-out tiles are refused** (`pipeline.check_allowed`) until `docs/prereg.md`
is marked `REGISTERED`.

## Pieces

| Module | What it does |
|---|---|
| `fetch.py` | Resumable, retrying downloads into `data/cache/` (connections that reset mid-transfer are resumed with HTTP Range) |
| `regions.py` | Tiles, 2 × 2-tile blocks, and the SHA-256-ordered calibration / held-out split (`scripts/make_split.py` prints it) |
| `prepare.py` | Copernicus GLO-30 mosaic and crop, EGM2008 → ellipsoidal heights (`h = H + N`, N read from the PROJ CDN grid), ESA WorldCover by majority resampling, UTM tile grids, ESRI ASCII grids with `.prj` |
| `pipeline.py` | `simulate_tile`: DEM crop with a 0.05° margin (outside terrain can shadow or lay over into the tile) → `bsar_render --geocode` with the real orbit and per-class materials → products on the tile's UTM grid. `compare`: every metric for one tile |
| `metrics.py` | Pure NumPy: IoU, 1-cell-tolerant F1, AUC (Mann–Whitney, ties averaged), 90 m block means, per-class Pearson r, phase-correlation offset |
| `evaluate.py` | `region-dem` (the ellipsoidal DEM for SNAP) and `run` (all tiles of a set vs. SNAP outputs, pooled as pre-registered, with pass/fail verdicts) |
| `dry_run.py` | End to end on calibration tiles with synthetic references whose answers are known |
| `materials_worldcover.txt` | WorldCover class → material (starting values; only `k` and `n` may be refitted, on calibration tiles) |
| `snap/s1_grd_beta0_tc.xml` | SNAP graph for the real-data side (orbit, noise removal, β⁰, terrain correction + layover/shadow mask) |

C++ support added for this pipeline: `bsar_eval passes` (closest approaches of
a real orbit to a point, with direction, look side and incidence), and
`bsar_render --geocode --classes --material-table` (range-Doppler geocoding
onto the DEM grid, with land-cover-dependent materials).

## Usage

```bash
pip install numpy rasterio
cmake --preset release && cmake --build --preset release

# Dry run on calibration tiles (no Sentinel-1 imagery needed)
python -m validation.dry_run --build build/release

# Real data (after SNAP has processed each scene with validation/snap/s1_grd_beta0_tc.xml)
python -m validation.evaluate region-dem --out build/validation/region_dem.tif
python -m validation.evaluate run --set calibration \
    --asc snap_asc.tif --asc-time 2024-08-14T17:15:51 \
    --desc snap_desc.tif --desc-time 2024-08-18T05:35:49 \
    --build build/release --out build/validation/calibration
```

## Status

| Piece | Status |
|---|---|
| DEM, geoid, WorldCover, UTM grids | Working; checked (geoid N = 50.5–52.7 m over the crop, matching the grid's own values) |
| Orbit download, pass finding | Working on real POEORB files; finds the S1A IW passes over the region (descending ~05:40, ascending ~17:20 UTC) |
| Simulation + geocoding per tile | Working, ~5–9 s per tile; masks tested against a closed-form ridge in `tests/test_geocode.cpp` |
| Metrics, pooling, SNAP-format loader | Working; unit tests in `validation/tests/`, plus identity and shift checks in the dry run |
| Dry run | Passes; report in [`docs/validation_dry_run.md`](../docs/validation_dry_run.md), and on every CI run |
| SNAP graph | **Written but never run.** To confirm with one real scene over calibration tiles: output band names (`Beta0_VV`, `layover_shadow_mask` expected) and mask encoding (bit 1 = layover, bit 2 = shadow expected); `evaluate run` takes both as options |
| Per-class `k`, `n` fitting | Not built yet; needs real calibration data |
| Real Sentinel-1 data | Needs a Copernicus Data Space account (not automated) |

## Dry-run quicklooks

Tile N46.60 E7.90 (Lauterbrunnen), simulated beta0 on the UTM grid with layover
tinted red and shadow blue. Opposite look directions put layover and shadow on
opposite slopes.

| Ascending, 2024-08-14 17:15 UTC (incidence ~33°) | Descending, 2024-08-18 05:35 UTC (~41°) |
|---|---|
| ![ascending](../docs/images/dryrun_asc_N46.60_E007.90.png) | ![descending](../docs/images/dryrun_desc_N46.60_E007.90.png) |
