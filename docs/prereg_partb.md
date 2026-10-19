# Pre-registration 3: a revised specular lobe (Part B)

> **Status: PROPOSAL, NOT YET REGISTERED.** A draft for review. To register:
> change anything you disagree with, fill in section 8, set this line to
> `REGISTERED`, and commit this file on its own; that commit's hash is the
> registration. After that, changes go only in "Amendments", dated and with a
> reason.
>
> **Nothing in the test region (section 3) may be rendered, downloaded for
> inspection, or looked at before registration.**

Part B of [`prereg_angular.md`](prereg_angular.md), registered separately as
that file's section 8 requires: once Part A had run there, its region no
longer counted as unseen. The revised model was developed on round-1 data only
([`partb_development.md`](partb_development.md)); round-2 data were not used
for development.

## 1. Questions

1. **Does the revised model fix the angular response on unseen data?** With
   the geometric-optics specular lobe scaled to 4.6% and the class parameters
   refitted on round 1, does the simulated response to local incidence match
   the real one on a new region, new scenes and new orbits?
2. **Does it keep what round 1 established?** Layover and shadow against
   SNAP, shadow against the real image, bare-rock correlation, geolocation.
3. **Does the current model's mismatch replicate a second time?** The
   current model is run on the same tiles as the reference.

## 2. Data

| Item | Choice | Why |
|---|---|---|
| Region | **Tyrol / Hohe Tauern, lon 11.0-13.0 E, lat 46.5-47.5 N** (Ötztal, Stubai and Zillertal Alps, Großglockner 3,798 m; the southern strip is South Tyrol) | Unseen (rounds 1 and 2 used 7-11 E, 46-47 N); strong relief; all four classes of interest; dense open weather data from GeoSphere Austria |
| Grid | UTM 32N at 30 m, as rounds 1 and 2, although 12-13 E lies in zone 33 | Keeps the registered SNAP graphs and tile code unchanged. At 13 E, 47 N the zone-32 scale error is about 0.07% (2 cm per 30 m cell) |
| Product, orbit, DEM, land cover, SNAP graphs | As rounds 1 and 2 (S1A IW GRDH VV, POEORB, Copernicus GLO-30 ellipsoidal, WorldCover 2021, ESA SNAP 14.0 with the round-1 graphs) | Only the region, scenes and model change |
| Scene selection rule | As rounds 1 and 2: datatake footprint ≥ 95% of the region, window 2024-07-15 to 2024-09-15, closest to 2024-08-15 that passes the snow check, one ascending and one descending acquisition | Same reasoning |
| Snow check | Same rule and thresholds as rounds 1 and 2 (wet hour ≥ 0.3 mm at any listed gauge; snowfall if a listed station at 2,500 m or higher is at or below +1.0 °C in a wet hour; reject on snowfall in the 72 h before, or rain at the acquisition hour). Data: GeoSphere Austria hourly station data (`klima-v2-1h`, open, no account): precipitation gauges at 1,500 m or higher in the region, temperature stations at 2,500 m or higher in or next to it, **listed in section 8 before any scene is selected** | MeteoSwiss has no stations here |

## 3. Tiles

0.1° tiles; **every one of the 200 tiles is a test tile.** Exclusions as round
1, section 3, with prereg 1 amendment 1 (a tile failing in either pass is
excluded from both). Development data for this test: round-1 data only.

## 4. Models, metrics and predictions

**Models (both frozen at registration):**

* **Current:** the round-1 registered model, `validation/materials_fitted.txt`
  (SHA-256 `d77d546c…`, LF; specular fraction 1).
* **Revised:** `validation/materials_partb_dev.txt` (SHA-256 `dbdf1747…`,
  LF): specular fraction 0.046 and k, n per class, fitted on round-1
  calibration tiles with `validation.fit --fit-specular-fraction`.

**Metrics:** A1-A4 and N exactly as `prereg_angular.md` section 4, computed by
`validation/angular.py` at the registered commit, for each model on the same
tiles and blocks.

| # | Metric | Revised model: threshold | Revised: prediction (round-1 development values) | Current model: prediction |
|---|---|---|---|---|
| A1 | Slope ratio, 20-60° | **0.75 to 1.33** per class | 0.80 to 1.15 (0.84-1.01) | 1.3 to 2.5 (as round 2) |
| A2 | Bias below 30° | **\|bias\| ≤ 2 dB** per class | within ±1.5 dB (−0.8 to +0.1) | ≥ +2 dB |
| A3 | Bias 40-60° | **\|bias\| ≤ 2 dB** per class | within ±1.5 dB (−0.3 to +0.3) | within ±2 dB |
| A4 | r, incidence removed | reported | reported; higher than the current model's for at least three of four classes | reported |
| N | Round-1 criteria G1-G3, R1, L1 | **each must pass** | pass | pass |

**Verdict.** Part B *succeeds* if, in both passes, the revised model meets the
A1, A2 and A3 thresholds for all four classes (tree cover, grassland, bare /
sparse vegetation, moss and lichen) and every N criterion passes. Anything else
is reported as a failure, class by class; no class is dropped after the fact.
The current model's numbers are reported next to the revised model's, and its
Part A pattern (A1 within 1.3-2.5 and A2 ≥ +2 dB for at least three classes in
both passes) is reported as a second replication. Every number is reported for
both models and both passes.

**Known risk, stated now:** bare rock's development A1 (0.84-0.86) is the
closest to its threshold, and in development rock alone preferred a larger
specular fraction (about 0.2). A single global fraction is kept on purpose.

## 5. What may be changed before registration, and on what

* Development data: round-1 data only. Nothing is fitted on the test region.
* The two models are fixed as above. Before registration the only code
  changes allowed are tooling: a GeoSphere reader for the snow check, the
  region in `validation/regions.py` and `select_scenes`, and the registration
  guard for this region (refused until this file is registered).

## 6. Procedure

1. **Before registration:** choose and list the snow-check stations (section
   8); add the tooling above with tests; run the scene selection (catalogue
   metadata only).
2. **Register:** fill in section 8 and commit this file on its own.
3. **Download and process** the scenes with the round-1 SNAP graphs; check the
   format on one block.
4. **Run once:** `validation.angular` for the current model and the revised
   model on all test tiles. Report every number in `docs/validation_partb.md`
   next to its threshold or prediction.

## 7. Known limitations, stated now

* One region, one season, two scenes.
* The Italian part (South Tyrol, south of the border at about 46.65-47.0 N)
  has no GeoSphere station; snow there is inferred from the Austrian
  stations, as round 2 did for its Italian part.
* The specular fraction is a single global number fitted on round-1 terrain;
  it does not say why the lobe is weaker than geometric optics predicts.
* k and n are not refitted on the test region, so A3 also tests how well
  round-1 parameters transfer.

## 8. Recorded at registration

* Snow-check stations: *to fill (before scene selection)*
* Scene IDs: *to fill*
* Pipeline commit: *to fill*

## Amendments

(none)
