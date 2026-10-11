"""Tiles, blocks and the geographic calibration / held-out split (docs/prereg.md, section 3).

Tiles are 0.1 x 0.1 degree cells named by their south-west corner, grouped into
0.2 x 0.2 degree blocks (2 x 2 tiles) so that neighbouring tiles, which share
mountains and therefore radiometry, never straddle the split. Blocks are ordered
by SHA-256 of "<seed>:<block>" (identical on every platform and Python version,
unlike random.shuffle) and alternately assigned, after forcing the blocks that
overlap areas already used during development into calibration.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass

# Region: Bernese and Valais Alps (Copernicus tiles N46E007 and N46E008).
LON_RANGE = (7.0, 9.0)
LAT_RANGE = (46.0, 47.0)
TILE = 0.1
BLOCK = 0.2
DEFAULT_SEED = 20261010

# Areas already looked at during development (must not be held out):
# the hillshade comparison crop (scripts/compare_hillshade.py).
DEV_AREAS = [(7.80, 46.50, 8.00, 46.62)]  # W S E N


@dataclass(frozen=True, order=True)
class Tile:
    """0.1-degree tile, identified by integer indices from the region's corner."""

    i: int  # longitude index
    j: int  # latitude index

    @property
    def bounds(self) -> tuple[float, float, float, float]:
        w = round(LON_RANGE[0] + self.i * TILE, 4)
        s = round(LAT_RANGE[0] + self.j * TILE, 4)
        return (w, s, round(w + TILE, 4), round(s + TILE, 4))

    @property
    def name(self) -> str:
        w, s, _, _ = self.bounds
        return f"N{s:05.2f}_E{w:06.2f}"

    @property
    def block(self) -> tuple[int, int]:
        return (self.i // 2, self.j // 2)


def blocks() -> list[tuple[int, int]]:
    nx = round((LON_RANGE[1] - LON_RANGE[0]) / BLOCK)
    ny = round((LAT_RANGE[1] - LAT_RANGE[0]) / BLOCK)
    return [(i, j) for j in range(ny) for i in range(nx)]


def block_bounds(b: tuple[int, int]) -> tuple[float, float, float, float]:
    w = LON_RANGE[0] + b[0] * BLOCK
    s = LAT_RANGE[0] + b[1] * BLOCK
    return (round(w, 4), round(s, 4), round(w + BLOCK, 4), round(s + BLOCK, 4))


def block_tiles(b: tuple[int, int]) -> list[Tile]:
    return [Tile(2 * b[0] + di, 2 * b[1] + dj) for dj in range(2) for di in range(2)]


def _overlaps(a: tuple[float, ...], b: tuple[float, ...]) -> bool:
    return a[0] < b[2] and b[0] < a[2] and a[1] < b[3] and b[1] < a[3]


def split(seed: int = DEFAULT_SEED) -> tuple[list[tuple[int, int]], list[tuple[int, int]]]:
    """(calibration blocks, held-out blocks), each sorted."""
    all_blocks = blocks()
    forced = [b for b in all_blocks if any(_overlaps(block_bounds(b), d) for d in DEV_AREAS)]
    rest = sorted((b for b in all_blocks if b not in forced),
                  key=lambda b: hashlib.sha256(f"{seed}:{b[0]},{b[1]}".encode()).hexdigest())
    calibration = list(forced)
    held_out: list[tuple[int, int]] = []
    for b in rest:  # forced blocks count towards calibration; keep the sizes equal
        (held_out if len(held_out) < len(calibration) else calibration).append(b)
    return sorted(calibration), sorted(held_out)


def calibration_tiles(seed: int = DEFAULT_SEED) -> list[Tile]:
    cal, _ = split(seed)
    return sorted(t for b in cal for t in block_tiles(b))


def held_out_tiles(seed: int = DEFAULT_SEED) -> list[Tile]:
    """Do not render, download for inspection, or look at these before registration."""
    _, held = split(seed)
    return sorted(t for b in held for t in block_tiles(b))
