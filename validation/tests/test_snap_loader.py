"""The SNAP reference loader: footprint and mask semantics found during the dry run."""

import numpy as np
import pytest

rasterio = pytest.importorskip("rasterio")

from validation import evaluate, prepare  # noqa: E402


def write_scene(path, beta0, mask, transform, crs):
    with rasterio.open(path, "w", driver="GTiff", width=beta0.shape[1], height=beta0.shape[0],
                       count=2, dtype="float32", crs=crs, transform=transform, nodata=0) as ds:
        ds.write(beta0.astype(np.float32), 1)
        ds.write(mask.astype(np.float32), 2)
        ds.set_band_description(1, "Beta0_VV")
        ds.set_band_description(2, "layover_shadow_mask")


def test_shadow_zeros_inside_swath_stay_valid(tmp_path):
    crs = prepare.UTM32N
    transform = rasterio.Affine(30.0, 0, 400000.0, 0, -30.0, 5160000.0)
    beta0 = np.full((20, 30), 0.1)
    beta0[:, :3] = 0.0           # outside the swath (left edge)
    beta0[8:12, 10:15] = 0.0     # radar shadow after thermal-noise removal
    beta0[5, 20] = -0.001        # over-subtracted noise
    mask = np.zeros_like(beta0)
    mask[8:12, 10:15] = 2        # shadow bit
    mask[2:4, 5:8] = 1           # layover bit
    path = tmp_path / "scene.tif"
    write_scene(path, beta0, mask, transform, crs)
    grid = prepare.Grid(transform, 30, 20, crs)
    ref = evaluate.load_snap_reference(path, grid, "Beta0_VV", "layover_shadow_mask", 1, 2)
    assert not ref["ref_valid"][:, :3].any()           # true no-data excluded
    assert ref["ref_valid"][8:12, 10:15].all()         # shadow zeros kept ...
    assert (ref["ref_beta0"][8:12, 10:15] == pytest.approx(1e-6))  # ... at the floor
    assert ref["ref_valid"][5, 20]
    assert ref["ref_shadow"][8:12, 10:15].all() and ref["ref_shadow"].sum() == 20
    assert ref["ref_layover"][2:4, 5:8].all() and ref["ref_layover"].sum() == 6
    assert ref["ref_valid"][0, 10]                     # mask value 0 is not no-data


def test_band_lookup_by_index_and_error(tmp_path):
    crs = prepare.UTM32N
    transform = rasterio.Affine(30.0, 0, 400000.0, 0, -30.0, 5160000.0)
    path = tmp_path / "scene.tif"
    write_scene(path, np.full((5, 5), 0.2), np.zeros((5, 5)), transform, crs)
    grid = prepare.Grid(transform, 5, 5, crs)
    ref = evaluate.load_snap_reference(path, grid, "1", "2", 1, 2)
    assert ref["ref_beta0"][2, 2] == pytest.approx(0.2)
    with pytest.raises(KeyError):
        evaluate.load_snap_reference(path, grid, "Sigma0_VV", "2", 1, 2)
