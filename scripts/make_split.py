#!/usr/bin/env python3
"""Print the deterministic geographic calibration / held-out split for
docs/prereg.md. The split itself lives in validation/regions.py.

Usage: python scripts/make_split.py [--seed 20261010] [--markdown]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from validation.regions import DEFAULT_SEED, block_bounds, block_tiles, split  # noqa: E402


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--seed", type=int, default=DEFAULT_SEED)
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
                print(" ", block_bounds(b), " ".join(t.name for t in block_tiles(b)))


if __name__ == "__main__":
    main()
