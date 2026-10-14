# Pre-registration 2: angular response of the scattering model

> **Status: REGISTERED on 2026-10-14 (Part A only).** The registration is the
> commit that set this line (this file changed on its own). From then on,
> changes go only in "Amendments", dated and with a reason; nothing above it is
> edited.
>
> **Nothing in the test region (section 3) may be rendered, downloaded for
> inspection, or looked at before registration.** Development uses round-1
> data only (section 5).

Follows [`prereg.md`](prereg.md) (round 1, registered at `0d5db40`) and the
exploratory analysis it led to ([`exploratory_incidence.md`](exploratory_incidence.md)),
which was done after the round-1 results were seen. This round turns that
exploration into a test on data nobody has looked at.

## 1. Questions

1. **Replication (Part A).** On a new region and new scenes, is the current
   model's angular response too steep, as the exploration suggested: does
   simulated β⁰ fall faster with local incidence than real β⁰, and is it too
   bright on slopes facing the sensor while matching at 40-60°?
2. **Fix (Part B).** Does a revised scattering model, developed on round-1
   data only, remove that mismatch on the new region without breaking what
   round 1 established (geometry, shadow, bare-rock correlation, geolocation)?

Part A can be registered and run on its own; Part B only once the revised
model exists (section 5). Part A is registered alone (section 8), so Part B's
metric definitions and thresholds are kept here for a later registration on a
further, unseen region; once Part A has run, this region no longer qualifies
as unseen for Part B.

## 2. Data

| Item | Choice | Why |
|---|---|---|
| Region | **Graubünden / Engadin, lon 9.0-11.0 E, lat 46.0-47.0 N** (Copernicus DEM tiles N46E009, N46E010) | Adjacent to round 1, so similar climate, land cover and sensor geometry, but no tile has been looked at. Strong relief (Bernina 4,049 m) and all four classes of interest (tree cover, grassland, bare rock, moss and lichen). Still UTM 32N |
| Product, orbit, DEM, land cover | As round 1 (S1A IW GRDH VV, POEORB, Copernicus GLO-30 ellipsoidal at 30 m, WorldCover 2021) | Only the region and the scenes change |
| Acquisitions | One ascending and one descending acquisition | Each slope seen from both sides |
| Scene selection rule | As round 1, section 2 (datatake footprint ≥ 95% of the region, window 2024-07-15 to 2024-09-15, closest to 2024-08-15 that passes the snow check), applied to this region with `validation.select_scenes`. Snow check: same rule and thresholds; stations: precipitation gauges at 1,500 m or higher in the region, temperature stations at 2,500 m or higher in or next to it, from MeteoSwiss open data, **listed in section 8 before any scene is selected** | Same reasoning as round 1; the station list has to change with the region |
| Round-1 scenes | Not used for any metric here | They are the development data |

## 3. Tiles

0.1° tiles as in round 1. **Every tile of the region is a test tile**: no part
of it is used for fitting or development, so no split is needed. Exclusions
as round 1, section 3, including amendment 1 (a tile failing in either pass is
excluded from both); the "mountainous" relief filter applies to the geometry
metrics only.

## 4. Metrics, thresholds and predictions

All on 90 m blocks (3 × 3 cells), per pass, pooled over test tiles.

**Blocks.** Class-pure (at least 7 of 9 cells of the class, as round 1) for
WorldCover 10 tree cover, 30 grassland, 60 bare / sparse vegetation and 100
moss and lichen; **no layover or shadow cell in the simulation or the SNAP
mask** (this differs from the exploration, which kept them; layover sums
several surfaces and would inflate the brightness of sensor-facing slopes in
both images). Snow and ice excluded.

**Local incidence** of a cell: angle between the DEM surface normal (central
differences on the 30 m UTM grid) and the direction to the sensor at the
tile centre's zero-Doppler time on the POEORB orbit, corrected for UTM grid
convergence (`validation/explore_incidence.py`, frozen at the registered
commit). A block's local incidence is the mean over its 9 cells.

**Uncertainty.** 95% intervals from 1,000 bootstrap resamples of tiles (seed
20261013). Pass/fail is decided on the point estimate; intervals are reported
next to it.

| # | Metric | Definition | Part A prediction (current model) | Part B threshold (revised model) |
|---|---|---|---|---|
| A1 | Angular slope ratio | Per class: slope of simulated dB against local incidence, divided by the real slope, by least squares on blocks with local incidence 20-60° | **1.3 to 2.5** for each of the four classes (exploration: 1.4-2.2) | **0.75 to 1.33** for each class |
| A2 | Sensor-facing bias | Per class: median of (simulated − real dB) over blocks with local incidence below 30° | **≥ +2 dB** for at least three of the four classes (exploration, 20-30° bin: +3.5 to +7.7 dB) | **\|bias\| ≤ 2 dB** for each class |
| A3 | Mid-range bias | Per class: median of (simulated − real dB) over blocks with local incidence 40-60° | **\|bias\| ≤ 2 dB** for each class (exploration: within about 1 dB; wider here because k and n were fitted on another region) | **\|bias\| ≤ 2 dB** for each class (set at 2 dB, not 1.5 dB, after the development run: tree cover sits at +1.4 to +1.5 dB on round-1 data, a level offset rather than an angular-shape error) |
| A4 | Correlation beyond incidence | Per class: Pearson r of simulated against real dB, each with its own local-incidence trend removed (mean in 2° bins of at least 30 blocks, interpolated; blocks outside the binned range dropped), as in the exploratory analysis | Reported, no prediction (exploration, with layover and shadow blocks kept: 0.44-0.63) | Reported; compared with the current model's value on the same blocks |
| N | Non-regression | Round-1 metrics G1, G2, G3, R1 and L1, computed exactly as round 1 on this region's test tiles | Reported for the current model (as a transfer check of round 1) | **Each must pass its round-1 threshold** |

**Part A verdict:** the exploration *replicates* if A1 holds for at least three
of the four classes in both passes and A2 holds in both passes. Anything else
is reported as a failure to replicate.

**Part B verdict:** the revised model *succeeds* if, in both passes, A1, A2 and
A3 meet their Part B thresholds for all four classes and every N metric passes.
Partial results are reported class by class; no class is dropped after the
fact.

Every number is reported for both models in Part B, and for both passes.

## 5. What may be developed, and on what

* **Development data:** round-1 region and scenes only (all its tiles, both
  calibration and held-out, are now considered seen).
* **Part A model:** frozen: the round-1 registered pipeline plus amendment 1,
  with `validation/materials_fitted.txt` (SHA-256
  `d77d546c4cde87c8378905d03527cbdef08607e03393db8fcf37d8347710e0f2`).
* **Part B model:** a change confined to the angular shape of
  `FacetScatteringModel` (near-normal / specular lobe, roughness handling).
  Permittivities, RMS slopes, geometry, geocoding and the metrics stay as they
  are. Its diffuse k and n are refitted with the round-1 procedure on the
  round-1 calibration tiles. The change, its commit and the refitted table are
  recorded in section 8 before registration.
* **Development evidence** required before registering Part B: the revised
  model's A1-A3 on round-1 data, labelled as development results (not tests),
  plus the round-1 dry run and the theory checks still passing.
* **Development baseline (current model, round-1 data, all 181 usable tiles):**
  [`angular_dev_round1.md`](angular_dev_round1.md), from `validation.angular`
  at `b31e8ba`. With these exact definitions the Part A criteria hold on
  round-1 data for all four classes in both passes: A1 1.38-2.33, A2 +3.2 to
  +7.8 dB, A3 within ±1.6 dB, A4 0.28-0.48. Bare rock and moss and lichen sit
  at A1 1.38-1.45, close to the predicted lower bound of 1.3; the bound was
  kept at 1.3 on purpose (the sharper claim). These are development results,
  not tests.
* **Never tuned on the test region:** anything.

## 6. Procedure

1. **Before registration:** list the snow-check stations (section 8); run the
   scene selection for the new region (catalogue metadata only); make the A1-A4
   metric code and test it on synthetic data and on round-1 data. For Part B,
   develop the revised model on round-1 data as in section 5.
2. **Register:** fill in section 8 (stations, scene IDs, commits, tables) and
   commit this file on its own.
3. **Download and process** the new scenes with the round-1 SNAP graphs; check
   format on one block, as round 1 did.
4. **Run once:** Part A (current model) and, if registered, Part B (revised
   model) on all test tiles. Report every number in `docs/validation_angular.md`
   next to its prediction or threshold.

## 7. Known limitations, stated now

* One region, one season, two scenes: the result describes late-summer Alpine
  terrain at the ellipsoid incidence of the selected scenes (recorded with the
  results), nothing wider.
* The per-tile look vector is taken at the tile centre (incidence varies by
  about 0.3° across a tile).
* WorldCover class purity is at 30 m; mixed pixels inside "pure" blocks remain.
* k and n are not refitted on the new region in Part A, so A3 also tests how
  well round-1 parameters transfer, not just the angular shape.

## 8. Recorded at registration

* Snow-check stations (chosen 2026-10-13 from the MeteoSwiss station
  metadata, before any scene selection for this region; every station has
  complete hourly data for 2024-07-12 to 2024-09-15):
  * Precipitation (wet hour if any reads ≥ 0.3 mm), every MeteoSwiss gauge at
    1,500 m or higher inside the region: San Bernardino (SBE, 1,639 m), Piz
    Martegnas (PMA, 2,668 m), Valbella (VAB, 1,568 m), Bivio (BIV, 1,856 m),
    Arosa (ARO, 1,878 m), Segl-Maria (SIA, 1,804 m), Weissfluhjoch (WFJ,
    2,691 m), Piz Corvatsch (COV, 3,294 m), Davos (DAV, 1,594 m), Samedan
    (SAM, 1,709 m), Passo del Bernina (BEH, 2,260 m), Buffalora (BUF, 1,971 m).
  * Temperature (snowfall if any is at or below +1.0 °C in a wet hour), every
    station at 2,500 m or higher inside the region: Piz Martegnas (PMA),
    Weissfluhjoch (WFJ), Piz Corvatsch (COV). Säntis (2,501 m) was left out:
    it is about 28 km north of the region, in a different massif (Alpstein).
  * Coverage gap: no MeteoSwiss station lies east of 10.3 E (the Austrian and
    Italian part of the region); snow there is inferred from the Swiss
    stations.
* Scene IDs (selection rule applied to the public catalogue on 2026-10-13 with
  `python -m validation.select_scenes --region graubuenden`, metadata only;
  `build/validation/round2/selected_scenes.json`):
  * Ascending, relative orbit 15, datatake 06C154, 99.3% of the region (two
    slices): `S1A_IW_GRDH_1SDV_20240826T171517_20240826T171542_055387_06C154_9E8E`,
    `S1A_IW_GRDH_1SDV_20240826T171542_20240826T171607_055387_06C154_04D0`.
    Zero-Doppler over 46.5 N 10.0 E: 2024-08-26 17:15:48 UTC, incidence 42.1°.
  * Descending, relative orbit 168, datatake 06BA07, 100% of the region (two
    slices): `S1A_IW_GRDH_1SDV_20240813T052706_20240813T052731_055190_06BA07_8CBE`,
    `S1A_IW_GRDH_1SDV_20240813T052731_20240813T052756_055190_06BA07_13E4`.
    Zero-Doppler over 46.5 N 10.0 E: 2024-08-13 05:27:35 UTC, incidence 41.5°.
  * Snow check: the ascending acquisition closest to 2024-08-15 (2024-08-14,
    06BAD3) fails, precipitating at the acquisition hour (no snowfall
    evidence). The selected ones pass with no snowfall evidence: ascending 12
    wet hours in the preceding 72 h (largest gauge total 9.3 mm, Buffalora),
    descending 4 wet hours (6.1 mm, Buffalora).
  * Both passes use different relative orbits from round 1 (88 and 66), so
    the sensor geometry is new as well as the terrain.
* Part A pipeline commit: `d3290709eb6bc133631e6bee716f2820fb70e39c` (metrics in
  `validation/angular.py`, local incidence in `validation/explore_incidence.py`),
  with `validation/materials_fitted.txt` (SHA-256 `d77d546c…`, LF) and ESA SNAP
  14.0 with the round-1 graphs. Every Part A number is produced with this
  commit, once.
* Part B: **not registered.** Running Part A means this region's images and
  results will have been seen, which could inform the development of a revised
  model. If Part B is pursued, it is therefore tested on a new region and new
  scenes under its own registration, not on this region.

## Amendments

1. **2026-10-14: region DEM for SNAP (tooling only).** At the registered
   commit, `python -m validation.evaluate region-dem` could only write the
   round-1 region's DEM, which SNAP's terrain correction needs for this region
   too. Added `--region` (default `bernese`, unchanged: regenerating the
   round-1 DEM with it is byte-identical). The round-2 DEM
   (`--region graubuenden`, 8.95-11.05 E, 45.95-47.05 N, 0.01% no-data) was
   made with it before any SNAP output or render of this region existed. No
   metric, render or threshold is affected. Found while preparing step 3.
   Also recorded here: the four downloaded products match the catalogue sizes,
   SHA-256 `2b9a7c39…` (06C154_9E8E), `1448e930…` (06C154_04D0), `a8d37ffd…`
   (06BA07_8CBE), `62d62a78…` (06BA07_13E4) (full values in the commit
   message).
