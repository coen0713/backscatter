# Round 3, Part B: a revised specular lobe on unseen terrain

Pre-registered in [`prereg_partb.md`](prereg_partb.md) (registered at
`961c5f0`, no amendments). Every number from the single registered run of
each model is reported here. Raw generated reports:
[current model](validation_partb_current_report.md),
[revised model](validation_partb_revised_report.md).

* Region: Tyrol / Hohe Tauern, 11.0-13.0 E, 46.5-47.5 N; all 200 tiles are test
  tiles, none looked at before registration. 165 kept in each pass; 35
  excluded under section 3, all for lying partly outside a footprint (20
  outside the ascending one, 15 outside the descending one; the region's
  west and east edges).
* Scenes: ascending 2024-08-09 17:07 UTC (relative orbit 117, 06B839),
  descending 2024-08-08 05:19 UTC (relative orbit 95, 06B755), two slices
  each, sizes matching the catalogue; ESA SNAP 14.0 with the round-1
  two-slice graph (Java heap capped at 8 GB; the graph is unchanged). Format
  checked on one block (11.8-12.0 E, 47.0-47.2 N) before the runs. Ellipsoid
  incidence at the tile centres: 37.6-46.6°.
* Models, both frozen at registration and run once each at `961c5f0`:
  **current** `validation/materials_fitted.txt` (round-1 fit, specular
  fraction 1) and **revised** `validation/materials_partb_dev.txt` (specular
  fraction 0.046, k and n refitted on round-1 calibration tiles).
* Product SHA-256: `51df7c1b…` (06B839_4E16), `d2c1b149…` (06B839_611A),
  `7279f7dc…` (06B755_6AB2), `501094e8…` (06B755_2B09) (full values in the
  commit message).

## Verdict

**Part B succeeds.** In both passes the revised model meets the A1, A2 and A3
thresholds for all four classes, and all five round-1 criteria (G1-G3, R1,
L1) pass. Every revised-model prediction held. The current model's mismatch
replicates a second time, as predicted.

## Revised model (ascending / descending; 95% tile-bootstrap intervals in the raw report)

| Class | A1 slope ratio (threshold 0.75-1.33; predicted 0.80-1.15) | A2 bias < 30° (threshold ±2 dB; predicted ±1.5) | A3 bias 40-60° (threshold ±2 dB; predicted ±1.5) |
|---|---|---|---|
| Tree cover | 1.00 / 0.97 | +0.03 / −0.07 dB | +0.02 / +0.07 dB |
| Grassland | 0.98 / 0.96 | +0.03 / −0.06 dB | +0.02 / −0.17 dB |
| Bare / sparse | 0.93 / 0.84 | −0.98 / −1.16 dB | −0.70 / −0.04 dB |
| Moss and lichen | 0.97 / 0.93 | −0.84 / −0.80 dB | −0.62 / −0.18 dB |

## Current model, same tiles and blocks (second replication of round 2)

| Class | A1 (predicted 1.3-2.5) | A2 (predicted ≥ +2 dB) | A3 (predicted ±2 dB) |
|---|---|---|---|
| Tree cover | 1.54 / 1.49 | +3.26 / +3.06 dB | +1.32 / +1.40 dB |
| Grassland | 2.28 / 2.22 | +8.00 / +7.92 dB | +0.15 / +0.02 dB |
| Bare / sparse | 1.57 / 1.38 | +3.32 / +3.01 dB | −1.10 / −0.53 dB |
| Moss and lichen | 1.44 / 1.35 | +3.19 / +3.24 dB | −0.99 / −0.49 dB |

The Part A pattern holds for all four classes in both passes. All five
round-1 criteria pass for the current model too.

## A4: correlation with the incidence trend removed (reported)

| Class | Current | Revised |
|---|---|---|
| Tree cover | 0.340 / 0.429 | 0.423 / 0.501 |
| Grassland | 0.488 / 0.415 | 0.577 / 0.556 |
| Bare / sparse | 0.432 / 0.333 | 0.440 / 0.425 |
| Moss and lichen | 0.511 / 0.320 | 0.551 / 0.435 |

The prediction (revised higher than current for at least three of four
classes) held for all four classes in both passes.

## Notes recorded before any result, and how they turned out

* **Bare rock was the stated risk** (development A1 0.84-0.86). It is again
  the weakest class: slightly too flat (A1 0.84 descending) and about 1 dB too
  dark facing the sensor, but inside every threshold. A single global
  specular fraction leaves rock slightly under-served, as development showed.
* **Wet spells before both scenes** (up to about 46 mm at one gauge in the
  72 h before): not tested separately. The current model's numbers here are
  close to its numbers in rounds 1 and 2, so no large effect is apparent.
* **A DEM error at Großglockner**, found while preparing the region DEM: the
  Copernicus GLO-30 summit is about 190 m too low over a few cells of tile
  N47.00_E012.60 (Großvenediger and Zuckerhütl check out within 10 m). Both
  simulation and SNAP use the same DEM, as pre-registered; the tile was kept.

## What this establishes

Across three pre-registered rounds (500 test tiles of Alpine terrain: 100 held-out in round 1, 200 in each of rounds 2 and 3):

* The geometry (layover, shadow, geolocation) agrees with ESA SNAP and the
  real images in every round.
* The original facet model's specular lobe is far too strong: its angular
  response is 1.35-2.5 times too steep on two unseen regions.
* Scaling the lobe to 4.6%, fitted on round-1 data alone, removes the
  mismatch on a third unseen region: slope ratios 0.84-1.00 and biases within
  1.2 dB, with nothing round 1 established lost.

It does not explain *why* the lobe is weaker than geometric optics predicts
(vegetation cover, sub-facet roughness and canopy attenuation are all
candidates), and a single global fraction is a description, not a physical
model. It covers late-summer Alpine terrain at 37-47° incidence only.
