"""Known-answer tests for the pre-registered metrics and the split."""

import numpy as np
import pytest

from validation import metrics, regions


def test_iou():
    a = np.array([[1, 1, 0], [0, 0, 0]], bool)
    b = np.array([[1, 0, 0], [1, 0, 0]], bool)
    valid = np.ones_like(a)
    assert metrics.iou(a, b, valid) == pytest.approx(1 / 3)
    valid[1, 0] = False  # ignored cells drop out of both masks
    assert metrics.iou(a, b, valid) == pytest.approx(1 / 2)
    assert np.isnan(metrics.iou(np.zeros_like(a), np.zeros_like(b), valid))


def test_tolerant_f1():
    a = np.zeros((20, 20), bool)
    a[5:15, 8] = True  # a thin line, like a shadow edge
    valid = np.ones_like(a)
    assert metrics.tolerant_f1(a, a, valid) == 1.0
    shifted = np.roll(a, 1, axis=1)
    assert metrics.iou(a, shifted, valid) == 0.0       # strict IoU collapses ...
    assert metrics.tolerant_f1(a, shifted, valid) == 1.0  # ... one cell is tolerated
    assert metrics.tolerant_f1(a, np.roll(a, 3, axis=1), valid) == 0.0
    assert np.isnan(metrics.tolerant_f1(np.zeros_like(a), np.zeros_like(a), valid))
    assert metrics.tolerant_f1(a, np.zeros_like(a), valid) == 0.0
    d = metrics.dilate(a, 1)
    assert d.sum() == 12 * 3 and d[4, 7] and not d[3, 8]


def test_auc():
    assert metrics.auc_lower_is_positive([1, 2], [3, 4]) == 1.0
    assert metrics.auc_lower_is_positive([3, 4], [1, 2]) == 0.0
    assert metrics.auc_lower_is_positive([1, 1], [1, 1]) == 0.5  # all ties
    # Against brute-force pair counting.
    rng = np.random.default_rng(0)
    pos, neg = rng.normal(0, 1, 200), rng.normal(0.7, 1, 300)
    pos[:20] = neg[:20]  # some exact ties
    pairs = (pos[:, None] < neg[None, :]).mean() + 0.5 * (pos[:, None] == neg[None, :]).mean()
    assert metrics.auc_lower_is_positive(pos, neg) == pytest.approx(pairs, abs=1e-12)


def test_block_reduce_and_fraction():
    v = np.arange(36, dtype=float).reshape(6, 6)
    valid = np.ones_like(v, bool)
    mean, ok = metrics.block_reduce(v, valid, k=3)
    assert mean.shape == (2, 2) and ok.all()
    assert mean[0, 0] == pytest.approx(v[:3, :3].mean())
    valid[0, 0] = False
    mean, ok = metrics.block_reduce(v, valid, k=3, min_fraction=1.0)
    assert not ok[0, 0] and ok[1, 1]
    assert mean[0, 0] == pytest.approx(v[:3, :3].ravel()[1:].mean())
    mask = np.zeros_like(v, bool)
    mask[:3, :3] = True
    assert metrics.block_fraction(mask, np.ones_like(mask), k=3)[0, 0] == pytest.approx(1.0)


def test_per_class_pearson():
    rng = np.random.default_rng(1)
    sim = rng.normal(size=(40, 40))
    classes = np.where(np.arange(40)[None, :] < 20, 10, 60) * np.ones((40, 1), int)
    real = np.where(classes == 60, sim + 0.1 * rng.normal(size=sim.shape), rng.normal(size=sim.shape))
    r = metrics.per_class_pearson(sim, real, classes, np.ones_like(sim, bool))
    assert r[60][0] > 0.95 and abs(r[10][0]) < 0.2
    assert r[60][1] == 800


@pytest.mark.parametrize("shift", [(0.0, 0.0), (2.0, -3.0), (0.4, 1.3)])
def test_phase_correlation_recovers_shift(shift):
    rng = np.random.default_rng(2)
    n = 128
    # A smooth random field, shifted exactly in Fourier space.
    field = np.real(np.fft.ifft2(np.fft.fft2(rng.normal(size=(n, n))) *
                                 np.exp(-0.002 * (np.add.outer(np.fft.fftfreq(n) ** 2,
                                                               np.fft.fftfreq(n) ** 2) * n * n))))
    ky = np.fft.fftfreq(n)[:, None]
    kx = np.fft.fftfreq(n)[None, :]
    moved = np.real(np.fft.ifft2(np.fft.fft2(field) * np.exp(-2j * np.pi * (ky * shift[0] + kx * shift[1]))))
    dy, dx = metrics.phase_correlation_offset(field, moved)
    assert dy == pytest.approx(shift[0], abs=0.25)
    assert dx == pytest.approx(shift[1], abs=0.25)


def test_split_is_the_registered_one():
    cal, held = regions.split(regions.DEFAULT_SEED)
    assert len(cal) == len(held) == 25
    assert not set(cal) & set(held)
    # Development area (hillshade crop) is never held out.
    dev_blocks = {t.block for t in regions.held_out_tiles()} & {(4, 2), (4, 3)}
    assert not dev_blocks
    # The exact held-out list quoted in docs/prereg.md (block south-west corners).
    quoted = ("46.4N 7.0E, 46.2N 7.2E, 46.0N 7.4E, 46.6N 7.4E, 46.8N 7.4E, 46.4N 7.6E, 46.6N 7.6E, "
              "46.0N 7.8E, 46.2N 7.8E, 46.8N 7.8E, 46.0N 8.0E, 46.6N 8.0E, 46.8N 8.0E, 46.0N 8.2E, "
              "46.8N 8.2E, 46.2N 8.4E, 46.6N 8.4E, 46.8N 8.4E, 46.2N 8.6E, 46.6N 8.6E, 46.8N 8.6E, "
              "46.0N 8.8E, 46.2N 8.8E, 46.4N 8.8E, 46.8N 8.8E")
    names = ", ".join(f"{regions.block_bounds(b)[1]:.1f}N {regions.block_bounds(b)[0]:.1f}E"
                      for b in held)
    assert names == quoted
    assert len(regions.calibration_tiles()) == 100
    assert regions.Tile(9, 5).name == "N46.50_E007.90"
