# Pre-registration: simulated vs. real Sentinel-1 (DRAFT, NOT YET REGISTERED)

> **Status: draft template.** Fill in every `TODO`, then commit this file on its
> own *before* rendering any real scene. The commit hash of that commit is the
> registration. Afterwards, any change goes in the "Amendments" section with a
> date and a reason; nothing above it is edited.

## Question

Does the geometric mode reproduce the layover and shadow geometry of real
Sentinel-1 IW acquisitions over mountainous terrain, and how well does the default
scattering model track real backscatter per land-cover class?

## Data

* Sensor and product: Sentinel-1 IW GRD, VV, POEORB orbits.
* Acquisitions: TODO (scene IDs, dates, ascending/descending).
* DEM: Copernicus GLO-30, converted to WGS84 ellipsoidal heights.
* Land cover: ESA WorldCover 2021.
* Reference masks: ESA SNAP layover/shadow mask (Terrain-Flattening /
  SAR-Simulation operator), same DEM and orbit. Software version: TODO.

## Tiles and split

* Tile size: TODO (e.g. 0.1° × 0.1°).
* Candidate region: TODO.
* Split rule (fixed before rendering anything): TODO. A deterministic rule
  such as "tile index parity", or a seeded shuffle with the seed recorded here.
* Calibration tiles: TODO (list).
* Held-out tiles: TODO (list).
* Exclusion criteria, decided now: TODO (e.g. > 20% DEM voids, snow cover in
  the acquisition month).

## Metrics and thresholds

| Metric | Computed on | Pass threshold | Rationale |
|---|---|---|---|
| Layover-mask IoU | held-out, mountainous tiles | TODO (e.g. ≥ 0.7) | TODO |
| Shadow-mask IoU | same | TODO | TODO |
| Pearson r of γ⁰ (dB), per WorldCover class | held-out tiles | TODO per class | TODO |
| Median geolocation offset of bright point targets | held-out tiles | TODO (e.g. < 1 pixel) | TODO |

Masks are compared in ground geometry after terrain correction, on the same
pixel grid. Per-class correlations require at least TODO pixels per class.

## What may be tuned, and on what

* May be tuned on calibration tiles only: material diffuse coefficients,
  exponents, RMS slopes.
* Never tuned: geometry code, orbit handling, thresholds above.

## Analysis plan

1. Render every calibration and held-out tile with the commit recorded at
   registration.
2. Fit the allowed parameters on calibration tiles only. Record the fitted values here.
3. Compute all metrics on held-out tiles once. Report all of them, including failures.
4. Report vegetation classes separately. They are expected to fail (no volume scattering).

## Amendments

(none)
