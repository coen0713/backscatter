"""The pre-registered metrics (docs/prereg.md, section 4), as pure NumPy.

All inputs are arrays on one common grid; `valid` masks select comparable cells.
"""

from __future__ import annotations

import numpy as np


def iou(a: np.ndarray, b: np.ndarray, valid: np.ndarray) -> float:
    """Intersection over union of two boolean masks over valid cells
    (G1, G2). Returns NaN if both masks are empty there."""
    a = a.astype(bool) & valid
    b = b.astype(bool) & valid
    union = np.count_nonzero(a | b)
    return float(np.count_nonzero(a & b) / union) if union else float("nan")


def dilate(mask: np.ndarray, radius: int = 1) -> np.ndarray:
    """Binary dilation with a (2r+1) x (2r+1) square."""
    out = mask.copy()
    h, w = mask.shape
    for dy in range(-radius, radius + 1):
        for dx in range(-radius, radius + 1):
            shifted = np.zeros_like(mask)
            shifted[max(dy, 0):h + min(dy, 0), max(dx, 0):w + min(dx, 0)] = \
                mask[max(-dy, 0):h + min(-dy, 0), max(-dx, 0):w + min(-dx, 0)]
            out |= shifted
    return out


def tolerant_f1(a: np.ndarray, b: np.ndarray, valid: np.ndarray, radius: int = 1) -> float:
    """Mask agreement that tolerates misregistration of up to `radius` cells
    (G1, G2 alongside strict IoU): precision = share of A within `radius` of B,
    recall = share of B within `radius` of A, F1 = their harmonic mean. Thin
    masks such as radar shadow lose most of their strict IoU to a one-pixel
    shift; this score does not. NaN if both masks are empty."""
    a = a.astype(bool) & valid
    b = b.astype(bool) & valid
    na, nb = np.count_nonzero(a), np.count_nonzero(b)
    if na == 0 and nb == 0:
        return float("nan")
    if na == 0 or nb == 0:
        return 0.0
    precision = np.count_nonzero(a & dilate(b, radius)) / na
    recall = np.count_nonzero(b & dilate(a, radius)) / nb
    return float(2 * precision * recall / (precision + recall)) if precision + recall else 0.0


def _average_ranks(x: np.ndarray) -> np.ndarray:
    order = np.argsort(x, kind="mergesort")
    xs = x[order]
    ranks = np.empty(len(x), dtype=np.float64)
    # Average ranks over runs of ties.
    boundaries = np.flatnonzero(np.diff(xs)) + 1
    starts = np.concatenate(([0], boundaries))
    ends = np.concatenate((boundaries, [len(x)]))
    for s, e in zip(starts, ends):
        ranks[order[s:e]] = 0.5 * (s + e - 1) + 1.0
    return ranks


def auc_lower_is_positive(positive: np.ndarray, negative: np.ndarray) -> float:
    """Area under the ROC curve for the rule "lower value => positive" (G3:
    real backscatter is lower where we predict shadow). 1 = perfect, 0.5 = chance.
    Mann-Whitney U with average ranks for ties."""
    positive = np.asarray(positive, dtype=np.float64).ravel()
    negative = np.asarray(negative, dtype=np.float64).ravel()
    n_p, n_n = len(positive), len(negative)
    if n_p == 0 or n_n == 0:
        return float("nan")
    ranks = _average_ranks(np.concatenate((positive, negative)))
    u_pos = ranks[:n_p].sum() - n_p * (n_p + 1) / 2.0  # P(positive > negative) * n_p n_n
    return float(1.0 - u_pos / (n_p * n_n))


def block_reduce(values: np.ndarray, valid: np.ndarray, k: int = 3,
                 min_fraction: float = 1.0) -> tuple[np.ndarray, np.ndarray]:
    """Mean over non-overlapping k x k blocks (30 m cells -> 90 m blocks for
    k = 3). A block is valid if at least `min_fraction` of its cells are."""
    h, w = (values.shape[0] // k) * k, (values.shape[1] // k) * k
    v = np.where(valid, values, 0.0)[:h, :w].reshape(h // k, k, w // k, k)
    m = valid[:h, :w].reshape(h // k, k, w // k, k)
    count = m.sum(axis=(1, 3))
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = v.sum(axis=(1, 3)) / count
    return mean, count >= min_fraction * k * k


def block_fraction(mask: np.ndarray, valid: np.ndarray, k: int = 3) -> np.ndarray:
    """Fraction of valid cells in each k x k block where `mask` is set."""
    frac, _ = block_reduce(mask.astype(np.float64), valid, k, min_fraction=0.0)
    return frac


def to_db(linear: np.ndarray, floor: float = 1e-6) -> np.ndarray:
    return 10.0 * np.log10(np.maximum(linear, floor))


def pearson(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=np.float64).ravel()
    b = np.asarray(b, dtype=np.float64).ravel()
    if len(a) < 3 or a.std() == 0 or b.std() == 0:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def per_class_pearson(sim_db: np.ndarray, real_db: np.ndarray, classes: np.ndarray,
                      valid: np.ndarray, min_count: int = 50) -> dict[int, tuple[float, int]]:
    """R1-R3: Pearson r of simulated vs. real dB per land-cover class, on
    blocks of a single (majority) class. Returns {class: (r, n_blocks)} for
    classes with at least `min_count` blocks."""
    out: dict[int, tuple[float, int]] = {}
    for c in np.unique(classes[valid]):
        m = valid & (classes == c)
        n = int(np.count_nonzero(m))
        if n >= min_count:
            out[int(c)] = (pearson(sim_db[m], real_db[m]), n)
    return out


def phase_correlation_offset(reference: np.ndarray, moving: np.ndarray) -> tuple[float, float]:
    """Sub-pixel (dy, dx) such that `moving` ~ `reference` shifted by (dy, dx)
    (L1). Both arrays must be the same shape and free of NaNs. Uses the
    normalised cross-power spectrum with a Hann window and a parabolic fit
    around the correlation peak."""
    if reference.shape != moving.shape:
        raise ValueError("phase_correlation_offset: shapes differ")
    win = np.outer(np.hanning(reference.shape[0]), np.hanning(reference.shape[1]))
    fa = np.fft.fft2((reference - reference.mean()) * win)
    fb = np.fft.fft2((moving - moving.mean()) * win)
    cross = fb * np.conj(fa)
    cross /= np.maximum(np.abs(cross), 1e-12)
    corr = np.fft.fftshift(np.real(np.fft.ifft2(cross)))
    py, px = np.unravel_index(np.argmax(corr), corr.shape)

    def refine(c_m: float, c_0: float, c_p: float) -> float:
        d = c_m - 2.0 * c_0 + c_p
        return 0.5 * (c_m - c_p) / d if d != 0 else 0.0

    h, w = corr.shape
    dy = refine(corr[(py - 1) % h, px], corr[py, px], corr[(py + 1) % h, px])
    dx = refine(corr[py, (px - 1) % w], corr[py, px], corr[py, (px + 1) % w])
    return float(py - h // 2 + dy), float(px - w // 2 + dx)


# ESA WorldCover 2021 class codes used by the prereg.
WORLDCOVER = {
    10: "tree cover", 20: "shrubland", 30: "grassland", 40: "cropland", 50: "built-up",
    60: "bare / sparse vegetation", 70: "snow and ice", 80: "permanent water bodies",
    90: "herbaceous wetland", 95: "mangroves", 100: "moss and lichen",
}
