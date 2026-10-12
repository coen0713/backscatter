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

## 2026-10-09: where the time goes (Phase 4)

**Backprojection was bound by `std::polar`, not memory.** Cache blocking
(pixel tiles x pulse blocks) gave 1.0-1.1x on one thread. The AVX2 kernel
does four pixels per step: double-precision range, an FMA-exact reduction of
the phase to [-pi/4, pi/4], Taylor sincos (error < 5e-12), and gathered
complex<float> interpolation with `addsub`. That gives 6-8x on one thread and
agrees with the scalar kernel to ~1e-7 relative. It is selected at runtime via
CPUID, and compiled per function with `target("avx2,fma")` rather than per
file, so no AVX2-encoded copy of a shared inline function can be picked by the
linker for non-AVX2 callers.

**Traversal is bound by memory latency, so wider nodes do not help.** The 8-wide
AVX2 BVH is within run-to-run noise of the 4-wide SSE one. Counting work per ray
on the 2.1M-triangle terrain explains it: 22 node visits and 2.2 triangle tests
per ray. At ~5.5 Mrays/s (coherent order) that is ~180 ns per ray, roughly 8 ns
per node visit, which is consistent with frequent cache misses into a 64 MB tree
rather than with arithmetic cost; random order (1.4 Mrays/s) is ~4x slower again.
(Not yet confirmed with hardware counters.) Larger SAH leaves were tried (intersection cost 0.25-1.0, up to
8 triangles per leaf) and were slower. The levers left are memory-side:
compressed or quantised nodes, treelet layouts, or ray packets/streams.

**Benchmarking on a hybrid CPU.** On the i9-13900H, single-thread results move
by about +-10% between runs depending on whether the thread lands on a P-core or
an E-core, and a run that overlapped other load showed 40% thread efficiency
instead of 88%. Benchmark numbers in the docs come from an idle machine on AC
power, and comparisons within ~10% are reported as ties.

## 2026-10-10: building the validation pipeline on calibration tiles

**Absolute orbit time was too coarse.** The pass finder found every closest
approach in a real POEORB file but Newton never converged. Orbit times were
seconds since 2000 (~7.7e8 s), where a double resolves only ~1.2e-7 s; Newton's
last steps (~3e-8 s) could not change `t`. Geometrically that is ~1 mm, but in
coherent mode 1 mm of range is ~0.2 rad of two-way C-band phase. `Orbit` now
works in orbit time (seconds since its first state vector, `epoch()`), with
absolute UTC only at the edges; Newton also stops when a step no longer changes
`t`. Synthetic test orbits started at t = 0, which is why the unit tests never
caught it; there is now a test with a 2024 epoch.

**A flaky network.** About one TLS connection in seven was reset on this
machine, for every client; GDAL's HTTP reader does not retry by default, so
remote COG reads failed or hung. Inputs are now downloaded once with a
resumable, retrying downloader into `data/cache/` and read locally.

**SNAP-reference loader, two traps found by the dry run's round trip.**
(1) Simulated deep shadow has beta0 = 0, and real data after thermal-noise
removal can be exactly 0 or negative there; treating 0 as no-data excluded every
tile through the footprint rule. The swath is now the span between the first and
last positive pixel of each row, with zeros inside kept at a -60 dB floor.
(2) The layover-shadow mask's value 0 means "neither"; using the band's no-data
value marked most of the scene invalid. Validity now comes from beta0 only.

**Strict IoU is the wrong score for thin masks.** Shifting a simulation against
itself by one 30 m cell (inside the 1-pixel geolocation tolerance) dropped
shadow IoU to 0.09-0.61 and layover IoU to 0.60-0.85. The proposed thresholds now
apply to an F1 that tolerates one cell; strict IoU is still reported. This was
changed before registration, from simulated data only.

**False alarms worth remembering.** A geoid check reported N = 44.6 m where
51.3 m was expected: the check compared two different DEM pixels on a cliff,
not the pipeline. A geocoding test failed because the "flat ground" probe sat in
passive layover (ground in front of a ridge shares range bins with its crest);
the test now asserts that passive layover explicitly.

**Materials matter even for a geometry test.** The first quicklooks showed lakes
as bright as meadows because the DEM got one material. The renderer now takes a
WorldCover class grid and a class -> material table.

## 2026-10-10: the per-class fit, and two versions that failed their own test

The fit is tested by parameter recovery: render a synthetic "real" image from a
truth table (three classes changed, 4-look speckle), fit from the starting
table, and require the truth back (|dn| <= 0.15, k within 10%).

**Version 1** predicted a block as `S + k F_n` from full two-bounce renders,
on blocks with 7 of 9 cells in the class. Forest and grassland came back, but
bare rock's k was 17% high, concrete 28% high, and water's n was 2.61 instead
of 2.0, although water and concrete were *unchanged* in the truth table. Two
bright forest cells in a shoreline block, or one cell of resampling bleed at a
class edge, dominate a dark class.

**Version 2** used interior blocks only (the class on all 9 cells and a 1-cell
ring). Contamination gone, but water still came back at n = 2.30. Rendering
the same experiment single-bounce recovered water exactly (k = 0.00300,
n = 2.00), and with double bounce 19% of the lake's signal was multi-bounce:
sensor -> water -> cliff -> sensor paths land at lake ranges but scale with the
cliff's k, which the model was crediting to water.

**Version 3** (current) attributes only single bounce to the block's class:
block = S1 + k D1(n) + M, with M (multi-bounce) rendered from the current table,
held fixed, and re-rendered after each of two fit rounds. All four classes with
enough blocks are recovered (forest 0.294/0.39 for 0.30/0.4, grassland
0.197/1.30 for 0.20/1.3, bare rock 0.0718/2.49 for 0.07/2.4, water 0.003/1.94
for 0.003/2.0).

Also changed: 90 m blocks now average linear power and then convert to dB
(multilooking) everywhere, instead of averaging dB values; renders are skipped
when the same command already completed, which makes long runs resumable.

## 2026-10-11: real scene selection, and a degeneracy that looked like a bias

**The draft selection rule could not be met.** Applied to the public catalogue,
no ascending Sentinel-1A acquisition in the window covers the whole region (best
97.2%, relative orbit 88; the east edge at 8.8-9.0 E is outside), and the
nearest descending one misses 0.04% of edge points. The rule now takes
acquisitions covering at least 95%, with uncovered tiles dropped by the existing
footprint exclusion; decided before registration and before any imagery was
obtained. An "acquisition" is a datatake: GRD products are 25 s slices, and the
descending scene is two of them. A useful cross-check fell out: the pass
finder's zero-Doppler times over the region centre, computed from the precise
orbits alone, lie inside the catalogue's slice time ranges for both scenes.

**Water stopped being recovered at the new geometry.** With the selected
ascending pass (42 deg instead of 33 deg), the fit returned water at k = 0.0035,
n = 2.51 instead of 0.003, 2.0. Without speckle it was exact on every tile; with
speckle it drifted. A lake is flat, so its blocks span only a few degrees of
incidence and k cos^n cannot separate k from n: the drifted pair predicts the
same backscatter as the truth to 0.06 dB. That is a degeneracy, not a bias, and
asking for the exact parameters was the wrong test. The fit now reports, per
class, how much its objective rises at n +- 0.5; below 0.1 dB n is "weakly
constrained" and is not interpreted, and the recovery check accepts a weakly
constrained class whose predictions match the truth.

## 2026-10-11: the selected scenes fail a snow check, so the fallbacks are used

The prereg rejects a scene with fresh snow above 2,500 m in the 3 days before it,
but left "the MeteoSwiss record" undefined. It is now a script
(`validation/snow_check.py`) on MeteoSwiss open hourly station data, fixed before
any imagery: an hour is wet if Guetsch or Pilatus records at least 0.3 mm, a wet
hour counts as snowfall if any station at or above 2,500 m is at or below +1 C,
and a scene fails on snowfall evidence in the preceding 72 h or rain at the
acquisition hour. (The SLF IMIS network measures snow height directly, but its
open API only reaches 7 days back.)

Both scenes closest to 15 August fail: a cold front on 17-18 August brought
precipitation with Jungfraujoch below freezing, and it was raining during the
18 August descending pass. The next closest eligible acquisitions pass and are
now selected: ascending 2024-08-07 (datatake 06B719, 97.3% coverage) and
descending 2024-08-06 (06B64D, two slices, full coverage). The geometry is the
same relative orbits as before, so the dry-run incidence angles do not change.

## 2026-10-12: SNAP's output format, confirmed on one calibration block

With the scenes downloaded (sizes match the catalogue), the SNAP graphs were run
for the first time, through `validation/snap_check.py`. It takes the production
graph and adds only a Subset after terrain correction, so a single calibration
block is computed and written; the subset is 0.01 deg inside the block, because
the UTM bounding box of a lat/lon rectangle pokes outside it and three of its
neighbours are held out.

What SNAP 14 actually writes: band 1 beta0, band 2 the mask, no band names and
no nodata value. The evaluation selected bands by name, so it now uses indices
(and the dry run's fake SNAP file mimics the real layout). The mask encoding was
checked without the simulator, so no pre-registered metric touched real data:
in both passes, mask value 1 sits on slopes rising about 35 deg along the look
direction (layover), value 2 on slopes falling about 43-45 deg (shadow), value 3
on both. SNAP fetched the EGM96 grid during the run, which raised the question
of whether it adds a geoid to our already-ellipsoidal DEM; re-running with
`externalDEMApplyEGM` true changes every pixel, so false is honoured.

The two-slice graph failed on first contact: SNAP refuses thermal-noise removal
after SliceAssembly. Noise is now removed per slice, then the slices are joined
and the orbit applied.
