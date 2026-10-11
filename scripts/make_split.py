#!/usr/bin/env python3
"""Deterministic geographic calibration / held-out split for docs/prereg.md.

Tiles are 0.1 x 0.1 degree cells named by their south-west corner. They are
grouped into 0.2 x 0.2 degree blocks (2 x 2 tiles) so neighbouring tiles, which
share mountains and therefore radiometry, never straddle the split. Blocks are
ordered by SHA-256 of "<seed>:<block>" (identical on every platform and Python
version, unlike random.shuffle) and alternately assigned, after forcing the
blocks that overlap areas already used during development into calibration.

Usage: python scripts/make_split.py [--seed 20261010] [--markdown]
"""

from __future__ import annotations

import argparse
import hashlib

# Region: Bernese and Valais Alps (Copernicus tiles N46E007 and N46E008).
LON_RANGE = (7.0, 9.0)
LAT_RANGE = (46.0, 47.0)
TILE = 0.1
BLOCK = 0.2

# Areas already looked at during development (must not be held out):
# the hillshade comparison crop (scripts/compare_hillshade.py).
DEV_AREAS = [(7.80, 46.50, 8.00, 46.62)]  # W S E N


def blocks() -> list[tuple[int, int]]:
    nx = round((LON_RANGE[1] - LON_RANGE[0]) / BLOCK)
    ny = round((LAT_RANGE[1] - LAT_RANGE[0]) / BLOCK)
    return [(i, j) for j in range(ny) for i in range(nx)]


def block_bounds(b: tuple[int, int]) -> tuple[float, float, float, float]:
    w = LON_RANGE[0] + b[0] * BLOCK
    s = LAT_RANGE[0] + b[1] * BLOCK
    return (round(w, 4), round(s, 4), round(w + BLOCK, 4), round(s + BLOCK, 4))


def overlaps(a: tuple[float, ...], b: tuple[float, ...]) -> bool:
    return a[0] < b[2] and b[0] < a[2] and a[1] < b[3] and b[1] < a[3]


def tile_names(b: tuple[int, int]) -> list[str]:
    w, s, _, _ = block_bounds(b)
    names = []
    for dj in range(2):
        for di in range(2):
            lat = s + dj * TILE
            lon = w + di * TILE
            names.append(f"N{lat:05.2f}_E{lon:06.2f}")
    return names


def split(seed: int) -> tuple[list[tuple[int, int]], list[tuple[int, int]]]:
    all_blocks = blocks()
    forced = [b for b in all_blocks if any(overlaps(block_bounds(b), d) for d in DEV_AREAS)]
    rest = sorted((b for b in all_blocks if b not in forced),
                  key=lambda b: hashlib.sha256(f"{seed}:{b[0]},{b[1]}".encode()).hexdigest())
    calibration = list(forced)
    held_out = []
    # Fill so both sets end up the same size (forced blocks count towards calibration).
    for b in rest:
        (held_out if len(held_out) < len(calibration) else calibration).append(b)
    return sorted(calibration), sorted(held_out)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--seed", type=int, default=20261010)
    p.add_argument("--markdown", action="store_true")
    args = p.parse_args()
    cal, held = split(args.seed)
    for name, group in (("Calibration", cal), ("Held-out", held)):
        if args.markdown:
            print(f"**{name}** ({len(group)} blocks, {4 * len(group)} tiles), block south-west corners:\n")
            print(", ".join(f"{block_bounds(b)[1]:.1f}N {block_bounds(b)[0]:.1f}E" for b in group) + "\n")
        else:
            print(f"{name}: {len(group)} blocks")
            for b in group:
                print(" ", block_bounds(b), " ".join(tile_names(b)))


if __name__ == "__main__":
    main()
