#!/usr/bin/env python3
"""Run `bsar_eval theory` and write docs/theory_results.md (quoted by the
README and docs/validation.md)."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

from provenance import ROOT, exe, header


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--build", type=Path, default=ROOT / "build" / "release")
    p.add_argument("--out", type=Path, default=ROOT / "docs" / "theory_results.md")
    args = p.parse_args()

    result = subprocess.run([str(exe(args.build, "apps", "bsar_eval")), "theory"],
                            capture_output=True, text=True)
    sys.stdout.write(result.stdout)
    sys.stderr.write(result.stderr)
    body = result.stdout.replace("## Theory checks\n\n", "", 1)
    args.out.write_text("# Theory check results\n\n" + header("scripts/run_theory.py", args.build) + body)
    print(f"wrote {args.out}")
    return result.returncode


if __name__ == "__main__":
    sys.exit(main())
