# Design notes

How `backscatter` turns a triangle scene and an orbit into a SAR image, and every
approximation it makes along the way.

## 1. Frames and precision

| Quantity | Frame | Type | Why |
|---|---|---|---|
| Scene geometry, BVH | local east-north-up (ENU), metres | `float` | memory bandwidth; SIMD box tests |
| Sensor position, path lengths, phase | same ENU frame | `double` | 4πR/λ at R ≈ 850 km is ~2·10⁸ rad; needs ~1 mm precision |
| Orbits | Earth-fixed (ECEF, from the EOF file) | `double` | mapped into ENU by `EnuFrame` |

A spaceborne sensor sits ~10⁶ m from the scene, where a `float` has ~6 cm
resolution. Rays from the sensor are therefore clipped to the scene bounding box
in double precision (`trace/sensor_rays.hpp`), and the BVH only ever sees a float
ray that starts just outside the geometry. Total path length is
`t_offset (double) + hit.t (float, small)`. The test "Sensor beyond float range
still bins correctly" checks the result against a double-precision range.

Earth-fixed velocities make zero-Doppler geometry correct without any Earth-rotation
term. Copernicus DEM heights are geoid (EGM2008) heights; `data/fetch_dem.py
--ellipsoidal` converts them, since orbits are ellipsoidal.

## 2. Acceleration structure

* **Binary BVH** (`geometry/bvh.cpp`). Binned SAH (16 bins by default; 2–64
  allowed) over triangle centroids, depth-first layout with 32-byte nodes, leaves of at most 4
  triangles. An object-median builder is kept as a baseline. Beyond depth 96 the
  builder switches to balanced splits, which bounds the depth (and so the fixed
  traversal stack) for any input.
* **Wide BVHs** (`geometry/wide_bvh.cpp`, `WideBvh<4>` and `WideBvh<8>`).
  Collapsed from the binary SAH tree by repeatedly opening the child with the
  largest surface area. Child boxes are stored lane-wise, so a single SIMD slab
  test checks all of them: SSE for 4-wide, AVX2 for 8-wide (chosen at run time
  from CPUID, with a scalar fallback). Hit children are pushed far-to-near so the
  nearest is popped first. Both widths perform about the same, because traversal
  is memory-latency bound (see `docs/devlog.md`); the scene uses the 4-wide tree,
  which needs nothing beyond SSE2.
* **Watertight intersection** (Woop, Benthin & Wald 2013). SAR looks at terrain at
  grazing angles, where ordinary Möller–Trumbore leaks through shared edges.
  The algorithm only works if `a*b - c*d` rounds both products separately: with FMA
  contraction (Clang's default when targeting a CPU with FMA), the two triangles
  sharing an edge stop computing exact negatives of each other's edge function, and
  rays slip through. The library is built with `-ffp-contract=off` as a public
  usage requirement. See `docs/devlog.md` for how this was found.

## 3. Geometric mode (`trace/geometric.cpp`)

For each azimuth line (time `t`, sensor `p`, velocity `v`):

1. **Rays in the zero-Doppler plane.** The plane through `p` perpendicular to `v`
   is spanned by `e1` (towards the scene centre) and `e2 = v̂ × e1`. The look
   angle sweeps the scene's bounding box in steps of
   `dθ = Δr / (rays_per_bin · R)`, so the ray spacing at the scene is a quarter of
   a range bin by default.
2. **Path tracing.** At each hit, with the normal facing the incoming ray:
   * bounce 1: `σ = σ₀(θ_local)` from `ScatteringModel::backscatter`;
   * bounce k > 1: if the sensor is visible (shadow ray), the bistatic coefficient
     `σ(n, incoming, to-sensor)` is weighted by the product of the Fresnel
     reflectances `|R|²` of the previous bounces;
   * then the ray reflects specularly. It stops after `max_bounces` or once the
     weight drops below 10⁻⁶.
3. **Radiometry.** A ray carries a beam tube of cross-section `A⊥ = R·dθ·Δaz`,
   which lands on area `A⊥ / cos θ_local`. The contribution `σ · A⊥ / cos θ_local`
   goes to the bin of the equivalent slant range `(path + |hit − p|) / 2`,
   normalised by the bin area `Δr·Δaz`. A flat plane then renders exactly
   `β₀ = σ₀ / sin θ` (tested to 5%).
4. **Anti-aliasing.** A direct return is spread uniformly over its footprint's
   slant-range extent `R·dθ·tan θ_local` instead of being binned as a point.
   Without this, the 4.77 rays per bin at 40° incidence alias into ±10% stripes.
   Multi-bounce paths stay point-binned: all paths in a dihedral have the same
   length, which is exactly why the corner line is so sharp.
5. **Layover and shadow masks.** For each line, neighbouring first hits on a
   continuous surface cover the slant-range interval between them. Two
   neighbouring hits count as continuous unless they are further apart than 8×
   the expected spacing `R·dθ / cos θ_local`; a larger gap is an occlusion edge. A sweep over interval
   endpoints gives the coverage multiplicity per bin:
   * **layover:** some part of the bin is covered more than once (> 5% of the bin);
   * **shadow:** more than half the bin, inside the swath, is covered by nothing.

Lines are independent, so threads never share output and the image is
bit-identical for any thread count.

**Approximations:** azimuth is binned at the line's zero-Doppler time for every
bounce order (exact for single bounce, approximate for paths that leave the
zero-Doppler plane); there is no antenna pattern or range spreading loss; and
polarisation is a scalar choice of Fresnel coefficient (no depolarisation).

## 4. Coherent mode (`trace/coherent.cpp`, `image/backprojection.cpp`)

1. **Scatterers.** Each facet gets `density × area` scatterers (stochastically
   rounded), placed uniformly with Philox4x32-10 keyed by
   `(seed, facet id, sample index)`. A visible scatterer carries
   `sqrt(σ₀ · area / n)`. No random phase is assigned: speckle comes from the
   random sub-wavelength positions alone.
2. **Multi-bounce paths.** From each scatterer, the specular reflection of the
   illuminating ray is traced on. Every further hit that sees the sensor becomes a
   path `(entry, exit, internal length)` with the geometric-mode weighting. Its
   range from sensor position `p` is `(|p − entry| + internal + |p − exit|) / 2`.
   Bounce points are found once (from the aperture centre) and held fixed. That
   is exact for corners parallel to the track and an approximation otherwise.
3. **Raw echoes.** For each pulse, every scatterer inside the azimuth beam
   (uniform: ±λ/2L; or the two-way sinc² main lobe) adds `a · e^{−i4πR/λ}` as a
   band-limited impulse at fractional delay `2R/c`. The impulse is a Kaiser-windowed
   sinc (32 taps, β = 6, 1024 tabulated phases). The impulse train is then convolved
   with the baseband LFM chirp by FFT. Pulses are independent (deterministic).
4. **Range compression.** Matched filtering by FFT correlation with the chirp,
   then 8× upsampling by zero-insertion in the middle of the spectrum, so
   backprojection can interpolate linearly. Normalised so a unit scatterer
   compresses to a unit peak.
5. **Backprojection.** For every pixel (any 3-D grid; the CLI drapes it on the
   scene surface), sum the compressed samples of all pulses at the pixel's range,
   phase-corrected by `e^{+i4πR/λ}`. Ranges and phases are in double precision.
   Kernels: pixel-major and cache-blocked (bit-identical to each other), and an
   AVX2/FMA kernel (four pixels per step, polynomial sincos after an FMA-exact
   phase reduction) that agrees with them to ~1e-7 relative and is used
   automatically when the CPU supports it. Each pixel always accumulates pulses
   in index order, so every kernel is bit-identical across thread counts.

**Approximations:** no range spreading loss (R⁻⁴) or elevation antenna pattern;
the stop-and-hop approximation (no motion during a pulse); visibility is
evaluated once, from the aperture centre.

## 5. Determinism

Every parallel loop writes disjoint outputs (one azimuth line, one pulse, one
pixel tile), and every reduction runs in a fixed order. Random numbers are a pure
function of `(seed, facet, index)`. `bsar_eval theory` checks that geometric and
coherent outputs are bit-identical on 1 thread vs. all hardware threads.
