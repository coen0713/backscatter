# Development log

Bugs worth remembering, written down while they were fresh.

## 2026-10-08: FMA contraction broke the watertight intersector

**Symptom.** The watertightness test fires rays from above a bumpy heightfield to
below it, aimed exactly at shared vertices and edge midpoints. 126 of 7,856 rays
passed straight through the surface. Brute force over every triangle missed them
too, so the BVH was not at fault.

**Diagnosis.** Dumping the edge functions for one failing near-vertical ray
through a cell diagonal showed the shared edge evaluated as
`V = −1.36036e-08` in one triangle and `W = −1.35857e-08` in the other. Woop et
al.'s argument needs these to be exact negatives (one triangle sees the ray just
inside, the other just outside). Both were negative, so both triangles rejected
the ray.

**Cause.** Zig's Clang targets the native CPU (which has FMA), and Clang contracts
`a*b − c*d` into `fma(a, b, −c*d)` by default. That rounds the two products
differently depending on operand order, which destroys the antisymmetry the
algorithm relies on.

**Fix.** `-ffp-contract=off` as a PUBLIC compile option of the library (the
intersector is inline in a header, so every consumer must agree). Zero leaks
afterwards. This is the grazing-angle seam leak the project plan predicted. The
lesson: "watertight" is a property of the arithmetic, not only of the algorithm.

**Also learned.** The first version of the test was wrong too. It aimed rays at
cell centres using bilinear interpolation, which is not on the triangulated
surface, and it counted rays that only graze a convex vertex as leaks. A ray
that touches a convex corner from above legitimately misses. The final test only
uses segments that cross the full height range of the surface, so any miss is a
real leak.

## 2026-10-08: rays-per-bin aliasing in the radiometry

**Symptom.** The flat-plane test expected β⁰ = 0.0730 and measured 0.0613.

**Diagnosis.** The row profile alternated between 0.0766 and 0.0613: bins were
receiving 5 or 4 rays. At 40° incidence there are 4/tan(40°) = 4.77 rays per bin,
so point binning aliases. The mean was right; individual bins were not.

**Fix.** Spread each direct return over its footprint's slant-range extent
`R·dθ·tan θ_local`, which integrates a flat surface exactly. Multi-bounce returns
stay point-binned on purpose: every path in a dihedral has the same length.

## 2026-10-08: chirp length off by one

`ceil(10e-6 * 40e6)` is 401, because the product is 400.00000000000006. Sample
counts derived from physical durations are now rounded to 1e-6 of a sample
before `ceil`.

## 2026-10-08: coherent mode was missing double bounce

The first coherent render of a building showed correct wall layover and shadow,
but no bright line at the wall foot: coherent scatterers were single-bounce only.
Scatterers now carry an entry point, an exit point and an internal path length, so
multi-bounce paths synthesise with range `(|p − entry| + L + |p − exit|) / 2`.
A unit test checks that the coherent dihedral focuses within 1 m of the wall foot,
more than 15 dB above the open ground.
