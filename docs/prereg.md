# Pre-registration: simulated vs. real Sentinel-1

> **Status: PROPOSAL, NOT YET REGISTERED.** Everything below is a concrete
> proposal for review. To register: change anything you disagree with, set the
> status line to `REGISTERED`, fill in the two "recorded at registration"
> fields (scene IDs, pipeline commit), and commit this file on its own. That
> commit's hash is the registration. From then on, changes go only in
> "Amendments", with a date and a reason; nothing above it is edited.
>
> **Nothing in the held-out set may be rendered, downloaded for inspection, or
> looked at before registration.** Building the pipeline on calibration tiles
> first is allowed and expected (see section 6).

## 1. Questions

1. **Geometry (primary).** Do the layover and shadow masks from the geometric
   mode agree with an independent implementation (ESA SNAP) on the same DEM and
   orbit, and do they predict where the *real* image is dark (radar shadow)?
2. **Radiometry (secondary).** After fitting the free parameters of the default
   scattering model on calibration tiles, how well does simulated backscatter
   track real backscatter per land-cover class on held-out tiles?
3. **Geolocation (secondary).** Is the simulated image registered to the real
   one to within a pixel?

## 2. Data

| Item | Choice | Why |
|---|---|---|
| Region | Bernese and Valais Alps, lon 7.0–9.0 E, lat 46.0–47.0 N (Copernicus DEM tiles N46E007, N46E008) | Strong relief (≈ 200 m at Lake Maggiore to 4,545 m at the Dom) makes layover and shadow common and large; good Sentinel-1 coverage from both orbit directions; the hillshade comparison already exercises this DEM |
| Product | Sentinel-1A IW GRDH, VV polarisation | Standard, free, 10 m pixels; VV is the co-pol channel the model represents |
| Orbit | POEORB precise orbit for each scene | ~5 cm accuracy; removes orbit error as an explanation |
| Acquisitions | **Two** scenes: one ascending and one descending relative orbit, each covering the whole region | Opposite look directions put layover and shadow on opposite slopes, so each slope is tested in both roles |
| Scene selection rule | For each direction: an *acquisition* is one datatake (all its GRD slices; footprint = union of slice footprints). Among acquisitions between 2024-07-15 and 2024-09-15 whose footprint covers **at least 95%** of the region, take the one closest to 2024-08-15 that passes the snow check, else the next closest. Snow check (`validation/snow_check.py`, MeteoSwiss open hourly data, fixed before any imagery was obtained): an hour is wet if Guetsch (2286 m) or Pilatus (2105 m) records >= 0.3 mm; a wet hour is snowfall evidence if any station at or above 2,500 m (Jungfraujoch, Gornergrat, Diablerets, Eggishorn, Les Attelas, Titlis) is at or below +1.0 °C. Reject the acquisition if there is snowfall evidence in the 72 h before it, or if it is precipitating at the acquisition hour. Tiles outside a selected footprint are dropped by the exclusions in section 3. Implemented by `python -m validation.select_scenes` | Fixes the choice before looking at any image; late summer minimises snow cover, which changes radiometry (not geometry); fresh snow and wet surfaces change C-band backscatter by several dB. The draft required full coverage; applying it to the catalogue showed no ascending acquisition in the window does (best 97.2%), so the threshold was set to 95% before registration and before any imagery was obtained |
| Scene IDs | Section 7 | |
| DEM | Copernicus GLO-30, converted to WGS84 ellipsoidal heights (`data/fetch_dem.py --ellipsoidal`), warped to UTM 32N at 30 m | Same DEM for simulation and for the reference; ellipsoidal heights to match the orbits |
| Land cover | ESA WorldCover 2021, resampled to the 30 m grid by majority | Defines the per-class strata |
| Buildings | None (the region is mostly rural; urban radiometry is out of scope here) | Keeps the test about terrain |

## 3. Tiles and split

* **Tiles:** 0.1° × 0.1° cells (≈ 7.7 km × 11 km) named by their south-west corner.
  The region holds 200 tiles.
* **Split unit:** 0.2° × 0.2° blocks of 2 × 2 tiles, so adjacent tiles that share
  a mountain cannot fall on both sides of the split.
* **Split rule:** `python scripts/make_split.py --seed 20261010`. Blocks are ordered
  by SHA-256 of `"<seed>:<block>"` and assigned alternately, after forcing every
  block that overlaps an area already used during development into calibration
  (currently the hillshade comparison crop, lon 7.80–8.00, lat 46.50–46.62).
  Result: 25 blocks (100 tiles) each:

  **Calibration** (block south-west corners): 46.0N 7.0E, 46.2N 7.0E, 46.6N 7.0E,
  46.8N 7.0E, 46.0N 7.2E, 46.4N 7.2E, 46.6N 7.2E, 46.8N 7.2E, 46.2N 7.4E, 46.4N 7.4E,
  46.0N 7.6E, 46.2N 7.6E, 46.8N 7.6E, 46.4N 7.8E, 46.6N 7.8E, 46.2N 8.0E, 46.4N 8.0E,
  46.2N 8.2E, 46.4N 8.2E, 46.6N 8.2E, 46.0N 8.4E, 46.4N 8.4E, 46.0N 8.6E, 46.4N 8.6E,
  46.6N 8.8E

  **Held-out** (block south-west corners): 46.4N 7.0E, 46.2N 7.2E, 46.0N 7.4E,
  46.6N 7.4E, 46.8N 7.4E, 46.4N 7.6E, 46.6N 7.6E, 46.0N 7.8E, 46.2N 7.8E, 46.8N 7.8E,
  46.0N 8.0E, 46.6N 8.0E, 46.8N 8.0E, 46.0N 8.2E, 46.8N 8.2E, 46.2N 8.4E, 46.6N 8.4E,
  46.8N 8.4E, 46.2N 8.6E, 46.6N 8.6E, 46.8N 8.6E, 46.0N 8.8E, 46.2N 8.8E, 46.4N 8.8E,
  46.8N 8.8E

* **Exclusions (decided now, applied per tile before any comparison):**
  * more than 5% DEM no-data → excluded entirely;
  * more than 50% WorldCover "permanent water" → excluded entirely (parts of lakes
    Thun, Brienz and Lucerne, and the north end of Lake Maggiore);
  * any part outside either scene's valid-data footprint → excluded entirely;
  * for the geometry metrics only, tiles must be **mountainous**: relief
    (max − min elevation) ≥ 1,000 m. Flat tiles have no layover or shadow to
    test, so they would only inflate agreement;
  * pixels of WorldCover "snow and ice" are excluded from the radiometry metric
    only (glaciers are bright and wet-snow-dependent; geometry is unaffected).

## 4. Metrics, thresholds and predictions

All comparisons are made on the same UTM 32N 30 m grid, per scene (ascending and
descending reported separately), over held-out tiles only.

| # | Metric | How | Pass threshold | Rationale |
|---|---|---|---|---|
| G1 | Layover agreement vs. SNAP | Our layover mask geocoded to the grid vs. SNAP Range-Doppler Terrain Correction "layover-shadow mask" (layover bit), same DEM and orbit, pooled over mountainous held-out tiles. **Thresholded on the 1-cell-tolerant F1** (a cell counts as matched if the other mask has a cell within 30 m); strict IoU is reported alongside | F1 ≥ 0.60 | Two independent implementations of the same geometry. Tolerant F1 rather than IoU: in the dry run, a 1-cell (30 m) shift alone, which L1 accepts, dropped strict layover IoU to 0.60–0.85 and shadow IoU to 0.09–0.61, while F1 stayed ≥ 0.978 (`docs/validation_dry_run.md`). SNAP's mask definition (e.g. whether passive layover is included) may differ from ours; differences will be analysed, not tuned away |
| G2 | Shadow agreement vs. SNAP | Same, shadow bit | F1 ≥ 0.60 | As above |
| G3 | Shadow vs. reality | AUC of real β⁰ (dB, 90 m blocks) for separating our predicted-shadow blocks (≥ 50% of cells in shadow) from blocks with **no** predicted shadow cell (partly shadowed blocks are left out); lower β⁰ = shadow | ≥ 0.85 | Tests the masks against the actual image, not another model: radar shadow sits at the noise floor (IW NESZ ≈ −22 dB) |
| R1 | Radiometry, bare rock / sparse vegetation | Pearson r of simulated vs. real β⁰ in dB of the 90 m block mean of linear β⁰ (multilooking ≈ 81 looks, so speckle is negligible), blocks that are class-pure (≥ 7 of 9 cells of the class), WorldCover 60, pooled over held-out tiles | ≥ 0.50 | Terrain-dominated surface scattering, which the model is built for |
| R2 | Radiometry, grassland | Same, WorldCover 30 | reported | Prediction: r between 0.3 and 0.6 (surface scattering plus some volume) |
| R3 | Radiometry, tree cover | Same, WorldCover 10 | reported | **Predicted to fail:** r < 0.3, because the model has no volume scattering. Reported in its own section, never averaged into a headline |
| L1 | Geolocation | Sub-pixel offset of real vs. simulated β⁰ (dB, 30 m cells, Hann-windowed phase correlation with parabolic peak fit) per held-out tile; median over tiles, east and north separately | ≤ 1 pixel (30 m) each | Orbit and DEM errors are far below this; a larger offset means a timing or datum bug |

**Comparison quantity.** Both sides use β⁰ (radar brightness), not γ⁰. Terrain
flattening would remove exactly the terrain effects being simulated. Real side
(SNAP, `validation/snap/s1_grd_beta0_tc.xml`): Apply Orbit File → Thermal Noise
Removal → Calibration (β⁰) → Range-Doppler Terrain Correction to UTM 32N at
30 m with the same ellipsoidal DEM, no radiometric normalisation. Our side:
geometric mode (2 bounces) at 10 m slant range / 10 m azimuth spacing, with
per-class materials from `validation/materials_worldcover.txt`, geocoded by
solving zero-Doppler for every DEM sample and sampling the slant-range image,
then resampled to the tile's UTM 30 m grid.

**Mask conventions** (matching SNAP's ground-geometry mask): *shadow* = the
ground cell cannot see the sensor (occlusion ray); *layover* = the cell's
slant-range bin receives returns from more than one part of the surface
(active and passive layover).

**Reference footprint.** Thermal-noise removal leaves shadow at exactly zero
(or below), so the swath is not "β⁰ > 0": per image row, everything between
the first and last positive pixel is inside it; zero or negative values inside
are kept at a −60 dB floor. The mask band's value 0 means "neither", never
no-data.

## 5. What may be tuned, and on what

* **Tunable, on calibration tiles only:** for each WorldCover class, the diffuse
  coefficient `k` and exponent `n` of `FacetScatteringModel`. They are fitted
  once, by minimising the median absolute dB error over calibration-tile 90 m
  blocks, and the fitted values are recorded below before any held-out tile is
  rendered. Procedure (`python -m validation.fit`, details in its docstring):
  * Blocks: interior blocks only (the class on all 9 cells and on a 1-cell ring
    around them), all cells valid, no layover or shadow cell; snow and ice
    excluded; both passes and all calibration tiles pooled. Classes with fewer
    than 100 such blocks keep their starting values.
  * Model: block = S1 + k · D1(n) + M, where S1 and D1(n) are single-bounce
    renders (diffuse off; k = 1 at exponent n on a 0.5 grid from 0 to 3,
    log-interpolated between nodes) attributed to the block's own class, and M
    is the multi-bounce part rendered with the current table and held fixed.
    Two rounds: fit, re-render M with the fitted table, fit again.
  * Search: for each n (0.05 steps, then golden-section), k by a log-spaced scan
    and golden-section refinement.
  * Identifiability is reported per class: the objective's rise at n ± 0.5
    with k re-optimised. Below 0.1 dB, n is "weakly constrained" (e.g. a flat
    lake seen over a few degrees of incidence), only the combination k·cosⁿ is
    determined, and the fitted n is not interpreted physically. The fitted pair
    is used as is either way: it is still the best predictor over that range.
  * Evidence the procedure recovers known parameters: the dry run fits a
    synthetic two-bounce "real" image with 4-look speckle from a truth table
    and recovers forest, grassland and bare rock within |Δn| ≤ 0.15 and 10% in
    k; water's k and n are degenerate at the selected geometry, and its fitted
    pair reproduces the true backscatter to 0.06 dB (`docs/validation_dry_run.md`). Two earlier versions failed that
    check and were changed before registration: 7-of-9 block purity let bright
    neighbours bias dark classes, and attributing multi-bounce signal to the
    block's own class biased water (lake → cliff → sensor paths).
* **Never tuned:** permittivities and RMS slopes (literature values), geometry
  code, orbit handling, the geocoding, the thresholds above, the split.
* **Bug fixes after registration:** allowed only if the bug is demonstrated on
  calibration tiles or synthetic tests. Each one is listed in Amendments with the
  commit, and all held-out metrics are recomputed and reported both before and
  after the fix.

## 6. Procedure

1. **Before registration** (allowed now): build the missing pipeline pieces using
   calibration tiles only, and dry-run end to end on calibration tiles.
   *Status: built and dry-run* (`validation/`, see `validation/README.md`;
   results in `docs/validation_dry_run.md`). *SNAP format confirmed*
   (`python -m validation.snap_check`) on both selected acquisitions over one
   calibration block only (7.8-8.0 E, 46.6-46.8 N, subset 0.01° inside its
   edges): band 1 = beta0, band 2 = mask with bit 1 = layover and bit 2 =
   shadow, no band names, so bands are selected by index. The mask encoding was
   checked against terrain slope along the look direction alone, not against
   the simulator, so no metric of section 4 has been computed on real data.
   The two-slice graph had to remove thermal noise per slice before slice
   assembly (SNAP refuses noise removal after it). The slice join was then
   checked on calibration block 7.2-7.4 E, 46.8-47.0 N, which the seam crosses
   at 46.90 N: the joined output has data in every row, each slice alone covers
   only its side, and the joined beta0 matches each slice to -0.05 / -0.08 dB
   mean (masks agree on 99.8% of pixels).
   *Scenes selected* by `python -m validation.select_scenes` (public catalogue
   metadata only, no imagery; see section 7). Their zero-Doppler times over
   46.5 N 8.0 E from the precise orbits (`bsar_eval passes`) fall inside the
   catalogue's slice times, an independent consistency check: ascending
   2024-08-07 17:24:01 UTC (incidence 42.4°), descending 2024-08-06 05:35:49 UTC
   (41.3°).
2. **Register:** record the scene IDs (catalogue only) and the pipeline commit
   hash below, set the status to `REGISTERED`, and commit.
3. Fit the tunable parameters on calibration tiles; record them in this file
   (as an amendment-free addition in section 8) and commit.
4. Render and compute every held-out metric **once**, with the registered commit
   plus the fitted parameters. Report every number, including failures, in
   `docs/validation.md`, with the prediction from section 4 next to each result.

## 7. Recorded at registration

* Scene IDs (from the selection rule applied to the catalogue on 2026-10-11;
  `build/validation/selected_scenes.json`):
  * Ascending, relative orbit 88, datatake 06B719, 97.3% of the region:
    `S1A_IW_GRDH_1SDV_20240807T172348_20240807T172413_055110_06B719_AF33`
  * Descending, relative orbit 66, datatake 06B64D, 100% of the region (two
    slices, joined with `validation/snap/s1_grd_beta0_tc_2slices.xml`):
    `S1A_IW_GRDH_1SDV_20240806T053519_20240806T053544_055088_06B64D_2F40`,
    `S1A_IW_GRDH_1SDV_20240806T053544_20240806T053609_055088_06B64D_8ECE`
  * Snow check: the acquisitions closest to 2024-08-15 (ascending 2024-08-19,
    06BD78; descending 2024-08-18, 06BCAB) both fail. A cold front on 17-18
    August brought precipitation with Jungfraujoch below 0 °C (first evidence
    2024-08-18 04:00 UTC, -0.3 °C), and it was precipitating during the
    descending acquisition. The next closest eligible acquisitions, listed above, pass:
    ascending 1.7 mm (Guetsch) / 11.7 mm (Pilatus) in the preceding 72 h with no
    station at or below +1 °C during wet hours; descending dry (0 / 0 mm).
* Pipeline commit: *to fill*

## 8. Fitted parameters

*(Filled in at step 3, before any held-out rendering.)*

## Amendments

(none)
