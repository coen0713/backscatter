"""Smoke tests for the Python bindings: shapes, dtypes and physics sanity."""

import numpy as np
import pytest

import backscatter


def test_render_synthetic_dihedral():
    r = backscatter.render_synthetic("dihedral", range_spacing=1.0, azimuth_spacing=4.0)
    intensity = r["intensity"]
    assert intensity.dtype == np.float32 and intensity.ndim == 2
    assert r["layover"].shape == intensity.shape and r["layover"].dtype == np.uint8
    assert len(r["bounce"]) == 3
    # The double-bounce line outshines all single-bounce returns.
    assert r["bounce"][1].max() > r["bounce"][0].max()


def test_render_heightmap_flat_has_no_layover_or_shadow():
    flat = np.zeros((33, 33), np.float32)
    r = backscatter.render_heightmap(flat, cell_size=30.0, incidence=40.0, bounces=1)
    assert r["layover"].sum() == 0
    assert r["shadow"].sum() == 0
    assert r["intensity"].max() > 0


def test_hillshade_flat_plane():
    shade = backscatter.hillshade(np.zeros((17, 17), np.float32), 10.0, elevation=45.0)
    assert shade.shape == (16, 16)
    np.testing.assert_allclose(shade, np.sin(np.radians(45.0)), rtol=1e-5)


def test_simulate_slc_speckle():
    slc = backscatter.simulate_slc(np.zeros((41, 41), np.float32), 1.0, pixel_size=2.0, density=8)
    assert slc.dtype == np.complex64
    intensity = np.abs(slc[4:-4, 4:-4]) ** 2
    enl = intensity.mean() ** 2 / intensity.var()
    assert 0.6 < enl < 1.6


def test_point_target_matches_theory():
    r = backscatter.point_target_check()
    assert r["range_3db"] == pytest.approx(r["expected_range_3db"], rel=0.05)
    assert r["azimuth_3db"] == pytest.approx(r["expected_azimuth_3db"], rel=0.05)
    assert r["range_pslr_db"] == pytest.approx(-13.26, abs=1.0)
