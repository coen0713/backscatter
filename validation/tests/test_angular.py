"""Round-2 metrics (docs/prereg_angular.md): known answers on synthetic data."""

import math

import numpy as np
import pytest

rasterio = pytest.importorskip("rasterio")

from validation import angular, explore_incidence, pipeline, prepare, regions  # noqa: E402


def blocks(theta, sim, real):
    return {"theta": np.asarray(theta, float), "sim": np.asarray(sim, float),
            "real": np.asarray(real, float)}


def test_slope_ratio_and_biases():
    th = np.linspace(0, 70, 701)
    real = -2.0 * th / 10
    sim = -4.0 * th / 10 + 6.0           # twice as steep, +6 dB at 0 deg
    m = angular.class_metrics(blocks(th, sim, real))
    assert m["A1"] == pytest.approx(2.0)
    assert m["slope_real_db_per_10deg"] == pytest.approx(-2.0)
    low = th < 30
    assert m["A2"] == pytest.approx(np.median((sim - real)[low]))
    mid = (th >= 40) & (th <= 60)
    assert m["A3"] == pytest.approx(np.median((sim - real)[mid]))
    assert m["A2_blocks"] == int(low.sum())


def test_too_few_blocks_give_nan():
    th = np.full(10, 25.0)
    m = angular.class_metrics(blocks(th, th, th))
    assert np.isnan(m["A1"]) and np.isnan(m["A2"]) and np.isnan(m["A3"])


def test_detrended_correlation_removes_the_incidence_trend():
    rng = np.random.default_rng(1)
    th = rng.uniform(10, 70, 20000)
    trend_r, trend_s = -0.3 * th, -0.6 * th     # strong, and of different steepness
    x = rng.normal(0, 1, th.size)               # what the model gets right beyond incidence
    m = angular.class_metrics(blocks(th, trend_s + x, trend_r + x + rng.normal(0, 1, th.size)))
    assert m["A4"] == pytest.approx(1 / math.sqrt(2), abs=0.03)
    m0 = angular.class_metrics(blocks(th, trend_s + rng.normal(0, 1, th.size),
                                      trend_r + rng.normal(0, 1, th.size)))
    assert abs(m0["A4"]) < 0.05                 # the raw r here is about 0.98


def test_bootstrap_is_reproducible_and_brackets_the_estimate():
    rng = np.random.default_rng(2)
    parts = []
    for _ in range(40):
        th = rng.uniform(10, 70, 200)
        real = -2.0 * th / 10 + rng.normal(0, 1, 200)
        parts.append(blocks(th, -3.0 * th / 10 + rng.normal(0, 1, 200), real))
    a = angular.bootstrap(parts, n=200, seed=7)
    b = angular.bootstrap(parts, n=200, seed=7)
    assert a == b
    point = angular.class_metrics(angular.pool(parts))["A1"]
    assert a["A1"][0] < point < a["A1"][1]
    assert a["A1"][0] > 1.2 and a["A1"][1] < 1.8   # truth 1.5


def _products(shape, layover=None, shadow=None, cls=30, beta0=0.1):
    grid = prepare.Grid(rasterio.Affine(30.0, 0, 500000.0, 0, -30.0, 5150000.0),
                        shape[1], shape[0], prepare.UTM32N)
    lay = np.zeros(shape, bool) if layover is None else layover
    sh = np.zeros(shape, bool) if shadow is None else shadow
    return pipeline.TileProducts(regions.Tile(0, 0), grid, np.full(shape, beta0), lay, sh,
                                 np.ones(shape, bool), np.zeros(shape), np.full(shape, cls))


def test_blocks_with_layover_or_shadow_are_dropped():
    shape = (9, 9)                             # 3 x 3 blocks
    lay = np.zeros(shape, bool)
    lay[0, 0] = True                           # one cell of block (0, 0)
    sim = _products(shape, layover=lay)
    ref_shadow = np.zeros(shape, bool)
    ref_shadow[8, 8] = True                    # one cell of block (2, 2), SNAP only
    ref = {"ref_beta0": np.full(shape, 0.2), "ref_valid": np.ones(shape, bool),
           "ref_layover": np.zeros(shape, bool), "ref_shadow": ref_shadow}
    b = angular.tile_blocks(sim, ref, np.full(shape, 40.0))
    assert b["sim"].size == 7                   # 9 blocks minus the two masked ones
    assert np.allclose(b["real"] - b["sim"], 10 * np.log10(2))
    assert np.all(b["cls"] == 30) and np.allclose(b["theta"], 40.0)


def test_other_classes_and_impure_blocks_are_not_used():
    shape = (6, 6)
    cls = np.full(shape, 30)
    cls[0:3, 0:3] = 20                         # shrubland: not one of the four classes
    cls[3:6, 3:6] = 30
    cls[3, 3] = cls[3, 4] = cls[3, 5] = 60     # 6 of 9 grassland: impure (needs 7)
    sim = _products(shape)
    sim.worldcover[:] = cls
    ref = {"ref_beta0": np.full(shape, 0.1), "ref_valid": np.ones(shape, bool),
           "ref_layover": np.zeros(shape, bool), "ref_shadow": np.zeros(shape, bool)}
    b = angular.tile_blocks(sim, ref, np.full(shape, 30.0))
    assert b["sim"].size == 2


def test_local_incidence_on_planes():
    # On the zone's central meridian grid north is true north, so the answer is exact.
    shape = (20, 20)
    grid = prepare.Grid(rasterio.Affine(30.0, 0, 500000.0, 0, -30.0, 5150000.0),
                        shape[1], shape[0], prepare.UTM32N)
    flat = np.zeros(shape)
    th = explore_incidence.local_incidence(flat, grid, 40.0, 90.0, 46.5, 9.0)
    assert np.allclose(th, 40.0, atol=1e-6)
    # Looking east (azimuth 90): a slope rising eastwards by 15 deg faces the sensor.
    x = np.arange(shape[1]) * 30.0
    ramp = np.tile(np.tan(np.radians(15.0)) * x, (shape[0], 1))
    th = explore_incidence.local_incidence(ramp, grid, 40.0, 90.0, 46.5, 9.0)
    assert np.allclose(th[2:-2, 2:-2], 25.0, atol=1e-6)
    th = explore_incidence.local_incidence(-ramp, grid, 40.0, 90.0, 46.5, 9.0)
    assert np.allclose(th[2:-2, 2:-2], 55.0, atol=1e-6)


def _res(a1, a2, a3):
    return {c: {"A1": a1, "A2": a2, "A3": a3} for c in angular.CLASSES}


def test_verdicts():
    good_a = {"ascending": _res(1.8, 4.0, 0.5), "descending": _res(1.8, 4.0, 0.5)}
    assert angular.part_a(good_a)["replicates"]
    one_pass_flat = {"ascending": _res(1.8, 4.0, 0.5), "descending": _res(1.0, 4.0, 0.5)}
    assert not angular.part_a(one_pass_flat)["replicates"]
    three_of_four = {d: dict(_res(1.8, 4.0, 0.5)) for d in ("ascending", "descending")}
    three_of_four["ascending"][10] = {"A1": 3.0, "A2": 1.0, "A3": 0.5}
    assert angular.part_a(three_of_four)["replicates"]
    nan_res = {d: _res(np.nan, np.nan, np.nan) for d in ("ascending", "descending")}
    assert not angular.part_a(nan_res)["replicates"]

    fixed = {d: _res(1.1, -1.0, 0.8) for d in ("ascending", "descending")}
    ok = {d: {"G1": True, "L1": True} for d in fixed}
    assert angular.part_b(fixed, ok)["succeeds"]
    regress = {"ascending": {"G1": True, "L1": False}, "descending": {"G1": True, "L1": True}}
    assert not angular.part_b(fixed, regress)["succeeds"]
    assert not angular.part_b({d: _res(1.1, -1.0, 1.6) for d in fixed}, ok)["succeeds"]


def test_round2_tiles_are_guarded():
    t = regions.region_tiles("graubuenden")[0]
    assert t.bounds == (9.0, 46.0, 9.1, 46.1)
    assert len(regions.region_tiles("graubuenden")) == 200
    if not pipeline.prereg_registered(pipeline.PREREGS["graubuenden"]):
        with pytest.raises(PermissionError):
            pipeline.check_allowed(t)
    with pytest.raises(ValueError):
        pipeline.check_allowed(regions.Tile(20, 0))   # past the round-1 region's east edge
