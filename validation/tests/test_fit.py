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


def test_identifiability_diagnostic():
    # Broad incidence range: n is constrained.
    y, s, f = synthetic_blocks(0.15, 1.2)
    assert not fit.fit_class(y, s, f)["n_weak"]
    # One incidence angle for every block (a flat lake): only k cos^n is
    # determined, so n must be flagged as weakly constrained.
    rng = np.random.default_rng(3)
    area = rng.uniform(0.8, 1.6, 2000)
    f1 = {g: area * 0.75 ** g for g in fit.N_GRID}
    r = fit.fit_class(0.15 * area * 0.75 ** 1.2, np.zeros(2000), f1)
    assert r["n_weak"]
    assert r["k"] * 0.75 ** r["n"] == pytest.approx(0.15 * 0.75 ** 1.2, rel=1e-3)


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


def test_tables_without_specular_fraction_are_unchanged(tmp_path):
    # The registered fitted table must read and rewrite byte-identically (no directive at 1).
    src = fit.ROOT / "validation" / "materials_fitted.txt"
    entries = fit.read_table(src)
    assert all(e.specular_fraction == 1.0 for e in entries)
    out = tmp_path / "t.txt"
    header = [line[2:] for line in src.read_text().splitlines() if line.startswith("# ")][:-1]
    fit.write_table(out, entries, "\n".join(header))
    assert out.read_text().splitlines() == src.read_text().splitlines()


def test_specular_fraction_directive_round_trip(tmp_path):
    base = fit.read_table(fit.pipeline.MATERIAL_TABLE)
    out = tmp_path / "t.txt"
    fit.write_table(out, fit.with_specular_fraction(base, 0.05), "test")
    assert "specular_fraction 0.05" in out.read_text()
    assert all(e.specular_fraction == 0.05 for e in fit.read_table(out))
    mixed = [base[0], fit.replace(base[1], specular_fraction=0.5)]
    with pytest.raises(ValueError):
        fit.write_table(out, mixed, "test")


def _lobe(cos):
    tan2 = (1 - cos ** 2) / cos ** 2
    return 0.2 * np.exp(-tan2 / (2 * 0.0625)) / (2 * 0.0625 * cos ** 4)


def test_recovers_a_global_specular_fraction():
    # The lobe matters on only ~20% of blocks, so the median objective sees it
    # only through realistic block-to-block scatter (here 9-look speckle,
    # about 1.4 dB): without noise the median is ~0 for any alpha and alpha
    # is not identified at all.
    rng = np.random.default_rng(5)
    pooled = {}
    for code, (k, n) in {30: (0.06, 0.5), 60: (0.2, 1.4)}.items():
        cos = rng.uniform(0.3, 1.0, 1500)
        area = rng.uniform(0.8, 1.6, 1500)
        spec = area * _lobe(cos)                    # S1: the full lobe
        f = {g: area * cos ** g for g in fit.N_GRID}
        y = (0.08 * spec + k * area * cos ** n) * rng.gamma(9, 1 / 9, 1500)  # truth: 8% of the lobe
        pooled[code] = fit.ClassBlocks([y], [spec], {g: [v] for g, v in f.items()}, [spec])
    alpha, profile = fit.fit_specular_fraction(pooled)
    assert alpha == pytest.approx(0.08, abs=0.02)
    at = {p["alpha"]: p["objective_db"] for p in profile}
    assert at[1.0] - min(at.values()) > 0.3      # the full lobe is clearly rejected
    fitted, report = fit.fit_all(pooled, fit.read_table(fit.pipeline.MATERIAL_TABLE), alpha)
    assert report[60]["k"] == pytest.approx(0.2, rel=0.15)
    assert all(e.specular_fraction == alpha for e in fitted)
