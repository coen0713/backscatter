# Round 2, Part A: angular response on an unseen region

Pre-registered in [`prereg_angular.md`](prereg_angular.md) (registered at
`269d169`, Part A only; amendment 1). Every number from the single registered
run is reported here. Raw generated report:
[`validation_angular_report.md`](validation_angular_report.md).

* Region: Graubünden / Engadin, 9.0-11.0 E, 46.0-47.0 N; all 200 tiles are test
  tiles, none looked at before registration. 185 kept in each pass; 15
  excluded under section 3, all for lying partly outside a footprint: 13 on
  the region's east edge (10.8-11.0 E) outside the ascending footprint, and
  N46.80_E009.00 and N46.90_E009.00 outside the descending one.
* Scenes: ascending 2024-08-26 17:15 UTC (relative orbit 15, 06C154),
  descending 2024-08-13 05:27 UTC (relative orbit 168, 06BA07), both two
  slices; reference ESA SNAP 14.0 with the round-1 two-slice graph. The
  format was checked on one block (9.8-10.0 E, 46.4-46.6 N) before the run.
* Model: current, frozen: round-1 pipeline, fitted table
  `validation/materials_fitted.txt` (fitted on round 1, not refitted here).
  Run at commit `7106eb0` (registered pipeline `d329070` plus amendment 1,
  tooling only).
* **The run was interrupted once.** The first attempt (2026-10-14) stopped
  after rendering 86 of 200 ascending tiles; the process was killed without
  error (it had been started detached from a shell that later closed). Metrics
  are computed only after every tile has rendered, so no metric existed or was
  seen. The run was restarted unchanged at the same commit and reused the
  86 completed renders (the cache keys on the exact render command and
  material table; rendering is deterministic). It completed on 2026-10-15 in
  2 h 06 min.

## Verdict

**Part A replicates.** In both passes, A1 lies within the predicted 1.3-2.5
for at least three of the four classes (ascending 4/4, descending 3/4), and
A2 is at least +2 dB for all four classes. A3 lies within ±2 dB for every
class in both passes, as predicted. All five round-1 criteria (G1-G3, R1, L1)
also pass on this region with the round-1 parameters.

## Results (ascending / descending, 95% tile-bootstrap intervals in the raw report)

| Class | A1 slope ratio (predicted 1.3-2.5) | A2 bias < 30° (predicted ≥ +2 dB) | A3 bias 40-60° (predicted within ±2 dB) | A4 r, incidence removed (reported) |
|---|---|---|---|---|
| Tree cover | 1.58 / 1.54 | +3.51 / +3.31 dB | +1.48 / +1.47 dB | 0.386 / 0.451 |
| Grassland | 2.13 / **2.52** | +8.25 / +8.74 dB | +0.22 / +0.49 dB | 0.477 / 0.376 |
| Bare / sparse | 1.48 / 1.51 | +3.54 / +3.38 dB | −0.69 / −0.70 dB | 0.511 / 0.301 |
| Moss and lichen | 1.38 / 1.61 | +3.91 / +4.43 dB | −0.38 / −0.43 dB | 0.523 / 0.278 |

**Outside its prediction:** grassland's A1 in the descending pass, 2.52
(interval 2.46-2.60), just above the predicted upper bound of 2.5: the
simulated angular response is even steeper relative to the real one than
predicted. Every other A1, A2 and A3 value lies within its prediction.

Slopes, dB per 10° of local incidence between 20° and 60° (real / sim):
tree cover −1.79 / −2.82 and −1.81 / −2.78; grassland −2.15 / −4.58 and
−1.82 / −4.60; bare −2.39 / −3.53 and −2.25 / −3.41; moss −2.59 / −3.58 and
−2.25 / −3.62.

## Comparison with the development baseline (round 1, not a test)

The values on the new region are close to the round-1 development values
([`angular_dev_round1.md`](angular_dev_round1.md)): A1 1.38-2.52 here vs.
1.38-2.33 there, A2 +3.3 to +8.7 dB vs. +3.2 to +7.8 dB, tree cover's A3
+1.5 dB in both. The bias on sensor-facing slopes is slightly larger on the
new region for every class.

## What this establishes (and what it does not)

* On late-summer Alpine terrain seen at 36.7-46.4° ellipsoid incidence (tile centres, both passes), the
  current facet model's response to local incidence is 1.4 to 2.5 times
  steeper than the real one, for all four natural classes, and it is 3.3 to
  8.7 dB too bright on slopes facing the sensor while matching to within
  1.5 dB between 40° and 60°. This was predicted from round-1 data and held on
  terrain, scenes and orbits nobody had looked at.
* It does not identify the cause. The exploration's candidates (too strong a
  near-normal / specular lobe, roughness too low) remain hypotheses; testing a
  revised model is Part B, which needs its own registration on a further
  unseen region (prereg 2, section 8).
