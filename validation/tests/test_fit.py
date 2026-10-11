"""The fitting maths: parameter recovery on synthetic blocks, and table I/O."""

import numpy as np
import pytest

from validation import fit


def synthetic_blocks(k, n, count=2000, noise_looks=None, seed=0):
    rng = np.random.default_rng(seed)
    cos = rng.uniform(0.2, 1.0, count)      # local incidence cosine per block
    area = rng.uniform(0.8, 1.6, count)     # geometric factor (slope, foreshortening)
    s = rng.uniform(0.0, 0.002, count)      # small specular part
    f = {g: area * cos ** g for g in fit.N_GRID}
    y = s + k * area * cos ** n
    if noise_looks:
        y = y * rng.gamma(noise_looks, 1.0 / noise_looks, count)
    return y, s, f


@pytest.mark.parametrize("k,n", [(0.12, 1.5), (0.30, 0.4), (0.07, 2.37), (0.2, 0.0)])
def test_recovers_parameters_exactly_without_noise(k, n):
    y, s, f = synthetic_blocks(k, n)
    r = fit.fit_class(y, s, f)
    assert r["n"] == pytest.approx(n, abs=0.01)
    assert r["k"] == pytest.approx(k, rel=0.01)
    assert r["median_abs_db"] < 0.01


def test_recovers_parameters_with_multilook_speckle():
    y, s, f = synthetic_blocks(0.15, 1.2, count=5000, noise_looks=81)
    r = fit.fit_class(y, s, f)
    assert r["n"] == pytest.approx(1.2, abs=0.1)
    assert r["k"] == pytest.approx(0.15, rel=0.05)


def test_diffuse_interpolation_is_exact_for_single_angles():
    _, _, f = synthetic_blocks(0.1, 1.0, count=50)
    cos_ratio = f[1.0] / f[0.0]
    np.testing.assert_allclose(fit.diffuse_at(f, 1.3), f[0.0] * cos_ratio ** 1.3, rtol=1e-9)
    assert fit.diffuse_at(f, 5.0) is f[3.0]  # clamped to the grid


def test_table_round_trip(tmp_path):
    base = fit.read_table(fit.pipeline.MATERIAL_TABLE)
    assert {e.code for e in base} >= {10, 30, 60, 80}
    out = tmp_path / "t.txt"
    fit.write_table(out, fit.uniform(base, 1.0, 0.5), "test")
    again = fit.read_table(out)
    assert [(e.code, e.eps_re, e.rms_slope) for e in again] == \
        [(e.code, e.eps_re, e.rms_slope) for e in base]
    assert all(e.k == 1.0 and e.n == 0.5 for e in again)


def test_small_classes_keep_starting_values():
    base = fit.read_table(fit.pipeline.MATERIAL_TABLE)
    y, s, f = synthetic_blocks(0.5, 1.0, count=fit.MIN_BLOCKS - 1)
    blocks = fit.ClassBlocks([y], [s], {g: [v] for g, v in f.items()})
    fitted, report = fit.fit_all({30: blocks}, base)
    assert report[30]["kept_start"]
    assert [e for e in fitted if e.code == 30][0].k == [e for e in base if e.code == 30][0].k
