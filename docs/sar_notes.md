# SAR in one page

Background notes for the code (Phase 0 deliverable). Symbols match `sensor/` and
`trace/coherent.hpp`.

## Geometry

A side-looking radar flies along a track (the **azimuth** direction) and points
sideways and down. **Slant range** `R` is the line-of-sight distance to a target;
**ground range** is its horizontal projection. The **incidence angle** θ is
measured from the local vertical at the target.

Two distortions follow directly from imaging by range rather than by angle:

* **Foreshortening and layover.** Slopes facing the radar are compressed in
  slant range. If a slope is steeper than the incidence angle (α > θ), its top is
  *closer* than its foot and the slope folds over the terrain in front of it:
  layover.
* **Radar shadow.** Slopes facing away and steeper than 90° − θ are never
  illuminated. Their slant-range interval receives no echo.

A wall standing on flat ground forms a **dihedral**. Paths that hit the ground
and then the wall (or the wall and then the ground) all have the same length as
the round trip to the corner, so they pile up into a bright line at the foot of
the wall: **double bounce**.

## Range resolution: chirps and matched filtering

A short pulse of duration τ resolves `cτ/2`, but carries little energy. SAR
instead transmits a long **linear FM chirp** `s(t) = exp(iπKt²)` of bandwidth
`B = Kτ`. Correlating the echo with the transmitted chirp (**range compression**,
a matched filter) concentrates the energy into a sinc of width `1/B` in time:

    slant-range resolution  δr = c / (2B)      (-3 dB width 0.886 c/(2B), unweighted)
    first sidelobe          -13.26 dB          (PSLR of a sinc)

C-band Sentinel-1 IW uses B ≈ 50 MHz, so δr ≈ 3 m. The `theory_radar()` used in
the checks uses 100 MHz, so δr = 1.5 m.

## Azimuth resolution: the synthetic aperture and Doppler

A real antenna of length `L` has beamwidth ≈ λ/L. At 850 km that is about 4 km
for Sentinel-1, which is useless for imaging. But each target stays inside the beam while
the platform flies a distance `L_s = R·λ/L`, the **synthetic aperture**. Over
that span the two-way phase `φ(t) = −4πR(t)/λ` follows a quadratic history
(the target's **Doppler** frequency sweeps through zero at closest approach).
Coherently summing the echoes along that history focuses the target, and the
achievable resolution is independent of range:

    stripmap azimuth resolution  δa = L / 2      (-3 dB width 0.886 L/2 with uniform weighting)

**Zero Doppler** is the instant when the sensor velocity is perpendicular to the
line of sight, i.e. the closest approach. `solve_zero_doppler` finds it by Newton
iteration on `f(t) = v(t)·(p(t) − x)`. Images are indexed by (zero-Doppler
azimuth time, slant range).

## Focusing by backprojection

Backprojection does the coherent sum directly in the time domain: for each output
pixel `x` and each pulse at position `p_k`,

    I(x) = Σ_k  s_k( R_k(x) ) · exp(+i 4π R_k(x) / λ),   R_k(x) = |p_k − x|

where `s_k` is the range-compressed echo of pulse k, interpolated at range `R_k`.
The phase term exactly undoes the propagation phase for a scatterer at `x`, so
those contributions add in phase while everything else averages out. It is
O(pixels × pulses): slow, but exact for any trajectory and any output grid,
including a grid draped on a DEM.

## Speckle

A resolution cell contains many scatterers at random sub-wavelength offsets. Their
phasors add as a 2-D random walk, so the complex pixel value is circular Gaussian,
amplitude is Rayleigh, and single-look **intensity is exponential**. Its standard
deviation equals its mean, so the **equivalent number of looks**
`ENL = mean² / variance = 1`. Averaging N independent looks gives ENL ≈ N. This
is a property of the physics, not noise, and the coherent mode reproduces it
without any random phases (see `docs/theory_results.md`).

## Radiometric conventions

* `σ⁰` (sigma-nought): radar cross-section per unit *ground* area.
* `β⁰` (beta-nought): radar brightness per unit *slant-range* area,
  `β⁰ = σ⁰ / sin θ` on flat terrain. The geometric mode outputs this.
* `γ⁰ = σ⁰ / cos θ`: often used for terrain-flattened products.
