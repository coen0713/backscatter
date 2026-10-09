#!/usr/bin/env python3
"""Run bsar_bench and write docs/benchmarks.md with commit and machine details."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

from provenance import ROOT, exe, header


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--build", type=Path, default=ROOT / "build" / "native")
    p.add_argument("--out", type=Path, default=ROOT / "docs" / "benchmarks.md")
    p.add_argument("--quick", action="store_true", help="small problem sizes (smoke test)")
    args = p.parse_args()

    cmd = [str(exe(args.build, "bench", "bsar_bench"))] + (["--quick"] if args.quick else [])
    result = subprocess.run(cmd, capture_output=True, text=True)
    sys.stdout.write(result.stdout)
    sys.stderr.write(result.stderr)
    if result.returncode != 0:
        return result.returncode
    args.out.write_text("# Benchmarks\n\n" + header("scripts/run_benchmarks.py", args.build) +
                        result.stdout)
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
