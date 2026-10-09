#!/usr/bin/env python3
"""Render the README figures into docs/images/ with bsar_render."""

from __future__ import annotations

import argparse
import shutil
import struct
import subprocess
import sys
import tempfile
import zlib
from pathlib import Path

from provenance import ROOT, exe

# (output name, bsar_render arguments, rendered file to keep)
FIGURES = [
    ("city_overlay", ["--scene", "city", "--range-spacing", "1", "--azimuth-spacing", "2"], "_overlay.png"),
    ("city_bounce2", ["--scene", "city", "--range-spacing", "1", "--azimuth-spacing", "2"], "_bounce2.png"),
    ("city_bounce1", ["--scene", "city", "--range-spacing", "1", "--azimuth-spacing", "2"], "_bounce1.png"),
    ("mountains_overlay", ["--scene", "mountains", "--range-spacing", "15", "--azimuth-spacing", "15"],
     "_overlay.png"),
    ("mountains_hillshade", ["--scene", "mountains", "--mode", "hillshade"], "_hillshade.png"),
    ("ridge_overlay", ["--scene", "ridge", "--range-spacing", "5", "--azimuth-spacing", "10"], "_overlay.png"),
    ("dihedral_coherent", ["--scene", "dihedral", "--mode", "coherent", "--platform", "airborne",
                           "--density", "2", "--pixel", "0.5"], "_intensity.png"),
]


def recompress_png(path: Path) -> None:
    """bsar_render writes stored (uncompressed) PNGs to avoid a zlib dependency;
    deflate them properly before they go into the repository."""
    data = path.read_bytes()
    pos, chunks, idat = 8, [], b""
    while pos < len(data):
        (length,) = struct.unpack(">I", data[pos:pos + 4])
        kind = data[pos + 4:pos + 8]
        body = data[pos + 8:pos + 8 + length]
        pos += 12 + length
        if kind == b"IDAT":
            idat += body
        else:
            chunks.append((kind, body))
    out = bytearray(data[:8])
    for kind, body in chunks:
        if kind == b"IEND":
            body_idat = zlib.compress(zlib.decompress(idat), 9)
            out += struct.pack(">I", len(body_idat)) + b"IDAT" + body_idat
            out += struct.pack(">I", zlib.crc32(b"IDAT" + body_idat))
        out += struct.pack(">I", len(body)) + kind + body + struct.pack(">I", zlib.crc32(kind + body))
    path.write_bytes(bytes(out))


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--build", type=Path, default=ROOT / "build" / "release")
    p.add_argument("--out", type=Path, default=ROOT / "docs" / "images")
    args = p.parse_args()
    render = exe(args.build, "apps", "bsar_render")
    args.out.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        done: dict[tuple[str, ...], str] = {}
        for name, render_args, suffix in FIGURES:
            key = tuple(render_args)
            if key not in done:
                prefix = str(Path(tmp) / name)
                print("+ bsar_render", " ".join(render_args))
                subprocess.run([str(render), *render_args, "--out", prefix], check=True)
                done[key] = prefix
            shutil.copy(done[key] + suffix, args.out / f"{name}.png")
            recompress_png(args.out / f"{name}.png")
            print(f"  -> {args.out / (name + '.png')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
