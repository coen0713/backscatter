# Pre-registration 2: angular response of the scattering model

> **Status: PROPOSAL, NOT YET REGISTERED.** A draft for review. To register:
> change anything you disagree with, fill in section 8, set this line to
> `REGISTERED`, and commit this file on its own; that commit's hash is the
> registration. After that, changes go only in "Amendments", dated and with a
> reason.
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
model exists (section 5).

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
| A3 | Mid-range bias | Per class: median of (simulated − real dB) over blocks with local incidence 40-60° | **\|bias\| ≤ 2 dB** for each class (exploration: within about 1 dB; wider here because k and n were fitted on another region) | **\|bias\| ≤ 1.5 dB** for each class |
| A4 | Residual correlation | Per class: Pearson r of (simulated − G) against (real − G) in dB, where G is a geometry-only render by the same model with every class given its fitted bare-rock material; removes local incidence and foreshortening alike | Reported, no prediction | Reported; compared with the current model's value on the same blocks |
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

* Snow-check stations: *to fill (before scene selection)*
* Scene IDs: *to fill*
* Part A pipeline commit: *to fill*
* Part B model change, commit and refitted table (SHA-256): *to fill, or
  "Part B not registered"*

## Amendments

(none)
