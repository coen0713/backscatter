# Part B development: the specular lobe (round-1 data only, not a test)

Everything here was done on round-1 data, which is development data for Part B
([`prereg_angular.md`](prereg_angular.md), section 5). None of it is a test:
the specular fraction and the class parameters were fitted on the same tiles
the metrics below are computed on. The test is a separate pre-registration on
an unseen region.

## Diagnosis

Round 2 showed the current model's response to local incidence is 1.4-2.5
times too steep and 3-9 dB too bright on sensor-facing slopes
([`validation_angular.md`](validation_angular.md)). Scaling the
geometric-optics specular lobe by a factor α, with each class's k and n
refitted by the registered procedure, removes the mismatch: at α = 1 the
round-1 numbers are reproduced, around α = 0.1 every class meets every Part B
threshold, and at α = 0 the response turns slightly too flat. So the lobe is
real but about twenty times too strong.

## Model change (commit `373e082`)

* `Material::specular_fraction` scales the lobe in `backscatter()` and
  `bistatic()`; mirror reflection on multi-bounce paths is unchanged. A
  material table sets it with one `specular_fraction <f>` line; tables
  without it behave and render exactly as before (the rebuilt renderer
  reproduces cached renders bit for bit, and the dry run's metrics and
  parameter recovery are identical).
* `python -m validation.fit --fit-specular-fraction` fits one fraction for all
  classes on the registered objective (median absolute dB error over the fit's
  blocks, pooled over classes, each class's k and n refitted at every α),
  then k and n as before. Without the flag the registered procedure runs
  unchanged.
* The objective only sees the lobe through block-to-block scatter: the lobe
  affects about a fifth of the blocks, so on noise-free synthetic data the
  median is flat in α. With realistic scatter a known α = 0.08 is recovered as
  0.079-0.083 (`validation/tests/test_fit.py`).

## Fit (round-1 calibration tiles, 96 of 100)

Fitted α = **0.046** (0.043 after the first iteration). The minimum is
shallow: every α between about 0.02 and 0.07 lies within 0.003 dB of it, so
the data say "about 5%". Report: [`partb_dev_fit.md`](partb_dev_fit.md);
table: `validation/materials_partb_dev.txt` (SHA-256, LF:
`dbdf1747bbeb29d5664e57bb97e0814f813304feaadf382dbd1ed89fd62fd067`).

| α | Pooled median abs. error |
|---|---|
| 1 (current model) | 1.72 dB |
| 0.3 | 1.00 dB |
| 0.1 | 0.91 dB |
| 0.046 | 0.90 dB |
| 0 | 0.91 dB |

Per class (k and n only → k, n and α): tree cover 1.74 → 0.73 dB, grassland
1.80 → 1.05, bare / sparse 1.77 → 1.48, moss and lichen 1.39 → 1.25; flat
classes barely change (water 1.45 → 1.43, cropland 1.20 → 1.17, built-up
2.30 → 2.30). Tree cover's n leaves the edge of the search range (0 → 0.91).

## Angular metrics with the revised table (all 200 round-1 tiles, in-sample)

Report: [`partb_dev_angular.md`](partb_dev_angular.md). Ascending / descending:

| Class | A1 slope ratio | A2 bias < 30° | A3 bias 40-60° |
|---|---|---|---|
| Tree cover | 1.01 / 1.01 | +0.03 / +0.10 dB | +0.00 / +0.15 dB |
| Grassland | 0.99 / 1.00 | −0.14 / −0.01 dB | −0.25 / −0.07 dB |
| Bare / sparse | 0.84 / 0.86 | −0.80 / −0.74 dB | +0.29 / +0.25 dB |
| Moss and lichen | 0.95 / 0.96 | −0.35 / −0.27 dB | +0.16 / +0.26 dB |

Every Part B threshold is met for every class in both passes, and all five
round-1 criteria still pass. Bare rock is now slightly too flat (A1 0.84-0.86,
−0.8 dB facing the sensor), consistent with rock alone preferring a larger
fraction (about 0.2) in the development scan. Theory checks pass unchanged.
