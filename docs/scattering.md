# Scattering model

`FacetScatteringModel` (`scene/material.hpp`) is deliberately simple, documented,
and swappable: integrators only talk to the `ScatteringModel` interface, so a
better model (e.g. a simplified IEM) can replace it without touching them.

## Terms

For a facet with unit normal `n`, local incidence θ, and material parameters
(ε, k, n_d, s):

**Diffuse term** (empirical, Lambert-like):

    σ⁰_d(θ) = k · cosⁿᵈ(θ)                 monostatic
    σ_d     = k · (cos θ_i · cos θ_s)^(n_d/2)   bistatic

**Specular lobe:** Kirchhoff scattering in the geometric-optics (stationary phase)
limit, for a surface with RMS facet slope `s`:

    σ⁰_s = |R(0)|² · exp(−tan²β / (2s²)) / (2 s² cos⁴β)

where β is the angle between `n` and the half vector of the incoming and outgoing
directions (β = θ for backscatter), and `R(0)` is the normal-incidence Fresnel
coefficient. Bistatically, the Fresnel coefficient is evaluated at the angle
between the incoming direction and the half vector.

**Specular reflection** (for multi-bounce paths) carries power `|R(θ)|²` with the
full complex-permittivity Fresnel coefficients:

    r_H = (cos θ − √(ε − sin²θ)) / (cos θ + √(ε − sin²θ))
    r_V = (ε cos θ − √(ε − sin²θ)) / (ε cos θ + √(ε − sin²θ))

A perfect conductor has `r_H = −1`, `r_V = +1`.

## Materials (C band, ~5.4 GHz)

| Name | ε (relative) | k | n_d | s | Source / note |
|---|---|---|---|---|---|
| `dry_soil` | 3.5 − 0.2j | 0.08 | 2 | 0.15 | Hallikainen et al. 1985, ~5% volumetric moisture |
| `wet_soil` | 15 − 3j | 0.15 | 2 | 0.15 | Hallikainen et al. 1985, ~25% moisture |
| `water` | 68 − 21j | 0.003 | 2 | 0.05 | Debye model, fresh water at 20 °C |
| `concrete` | 5.5 − 0.6j | 0.06 | 1.5 | 0.10 | Ulaby & Long 2014 |
| `asphalt` | 4.0 − 0.3j | 0.02 | 2 | 0.08 | Ulaby & Long 2014 |
| `metal` | perfect conductor | 0.01 | 1.5 | 0.05 | |

The diffuse coefficients are chosen so that σ⁰ at 35° falls in typical C-band
ranges (dry soil ≈ −13 dB, calm water < −20 dB). The unit tests check these
values and the Fresnel limits (Brewster angle, grazing incidence, normal incidence).

## Known failure modes (where this model is wrong)

These are reported separately in the validation rather than averaged away:

* **Vegetation.** Forest and crops scatter from a volume, not a surface. The model
  has no volume term, so vegetated areas will be too dark and too sensitive to
  terrain slope.
* **Energy bookkeeping.** The diffuse and specular terms are not normalised
  against each other, and the coherent specular reflection is not reduced by
  surface roughness (no Rayleigh roughness factor). Multi-bounce returns from rough
  surfaces are therefore overestimated.
* **Polarimetry.** Only a scalar HH or VV Fresnel coefficient is used, with no
  depolarisation and no cross-pol channel.
* **Wet snow, urban clutter, and Bragg scattering from the ocean** are out of scope.

## References

* F. T. Ulaby and D. G. Long, *Microwave Radar and Radiometric Remote Sensing*,
  University of Michigan Press, 2014.
* M. T. Hallikainen et al., "Microwave dielectric behavior of wet soil, part I",
  *IEEE TGRS* 23(1), 1985.
* S. Auer et al., "Ray-tracing simulation techniques for understanding
  high-resolution SAR images", *IEEE TGRS* 48(3), 2010 (RaySAR).
