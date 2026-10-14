# Exploratory: does the model get more than terrain right?

**Not pre-registered.** This follows up the two failed predictions in
[`validation.md`](validation.md) (R2 grassland, R3 tree cover correlated far
better than predicted). It was designed after the held-out results were seen,
so nothing here is a test; it is a description that suggests what to
pre-register next. Generated tables:
[`exploratory_incidence_report.md`](exploratory_incidence_report.md); script:
`validation/explore_incidence.py`.

**Method.** Same 85 held-out tiles and the same class-pure 90 m blocks as the
registered R metrics (the raw r values reproduce them exactly). For every 30 m
cell, local incidence is computed from the DEM normal and the look vector at
the tile's zero-Doppler time (POEORB orbit). Per class and pass, the mean dB in
2° local-incidence bins is removed from the real and simulated block values
separately, and the residuals are correlated. The script reads the existing
held-out renders and never renders.

## Results (ascending / descending)

| Class | r, registered | r, incidence removed | Real variance explained by local incidence | Slope 20-60°, dB per 10°: real | sim |
|---|---|---|---|---|---|
| Tree cover | 0.871 / 0.917 | 0.443 / 0.628 | 78% / 79% | −2.0 / −2.1 | −3.1 / −3.2 |
| Grassland | 0.860 / 0.861 | 0.564 / 0.562 | 67% / 66% | −2.2 / −2.3 | −4.8 / −4.9 |
| Bare / sparse | 0.876 / 0.883 | 0.542 / 0.553 | 73% / 73% | −2.7 / −2.7 | −4.1 / −4.1 |
| Moss and lichen | 0.900 / 0.894 | 0.481 / 0.524 | 79% / 75% | −2.6 / −2.9 | −3.8 / −4.2 |
| Cropland | 0.277 / 0.318 | 0.146 / 0.174 | 7% / 10% | −1.6 / −2.0 | −4.5 / −3.3 |
| Built-up | 0.240 / 0.253 | 0.180 / 0.225 | 4% / 3% | −1.7 / −1.4 | −1.9 / −2.0 |
| Water | 0.707 / 0.717 | 0.620 / 0.747 | 22% / 7% | (narrow range) | |

Median dB by local incidence, grassland, ascending (real / sim): 0-10°
2.3 / 13.3, 10-20° −1.0 / 8.8, 20-30° −4.9 / 2.8, 30-40° −7.8 / −4.1,
40-50° −9.7 / −9.7, 50-60° −11.1 / −11.9, 60-70° −12.6 / −13.2. The other
terrain classes and the descending pass show the same pattern (full tables in
the report).

## What this suggests

1. **The registered R metrics mostly measured terrain.** Local incidence alone
   explains two thirds to four fifths of the real block variance for every
   vegetated and rocky class. That is why R2 and R3 came out far above their
   predictions: they were never sensitive to the scattering mechanism the
   predictions were about.
2. **The model does get more than local incidence right.** With the incidence
   trend removed, correlations of 0.44-0.63 remain. Part of that is other
   geometry (foreshortening of the ground area, multi-bounce near cliffs), so
   this is still not a pure test of the scattering model.
3. **Forest is not the outlier.** Its residual correlation (0.44 / 0.63) is in
   the same range as bare rock (0.54 / 0.55), and its real angular response is
   only slightly flatter than grassland's (−2.0 vs −2.2 dB per 10°). In these
   data, the missing volume term does not set forest apart.
4. **The model's angular response is too steep for every natural class.**
   Between 20° and 60° the simulation falls 1.4-2.2 times faster than reality,
   and on slopes facing the sensor (local incidence under 30°) it is 4-11 dB
   too bright, worst for grassland. Above 40° the two agree to about 1 dB. The
   fit could not fix this: its objective is a median over all blocks, most of
   which lie at 35-60°, and a single k·cosⁿ term cannot both match that range
   and stay flat towards normal incidence. Candidates: too much specular /
   near-normal lobe in `FacetScatteringModel`, or roughness too low.
5. **Cropland and built-up** are flat valley floors: the model puts most of its
   variance on incidence (59-79%), the real image almost none (3-10%), so what
   varies in reality there is not geometry and the model has nothing to match
   it with.

## What would test it properly

A pre-registered follow-up, on new scenes or a fresh held-out split: (a) the
angular-response slope per class, real vs. simulated, with a threshold; (b)
residual correlation after removing a geometry-only render (uniform material),
which also removes foreshortening; (c) the low-incidence bias (local incidence
< 30°) as its own metric, since it is where the model fails most.
