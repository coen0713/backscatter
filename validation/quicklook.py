"""Small PNG quicklooks without extra dependencies."""

from __future__ import annotations

import struct
import zlib
from pathlib import Path

import numpy as np


def write_png(path: Path, image: np.ndarray) -> None:
    """8-bit grayscale (H x W) or RGB (H x W x 3) PNG."""
    img = np.clip(image, 0, 255).astype(np.uint8)
    color_type = 2 if img.ndim == 3 else 0
    raw = b"".join(b"\x00" + row.tobytes() for row in img)

    def chunk(kind: bytes, body: bytes) -> bytes:
        return struct.pack(">I", len(body)) + kind + body + struct.pack(">I", zlib.crc32(kind + body))

    ihdr = struct.pack(">IIBBBBB", img.shape[1], img.shape[0], 8, color_type, 0, 0, 0)
    path.write_bytes(b"\x89PNG\r\n\x1a\n" + chunk(b"IHDR", ihdr) +
                     chunk(b"IDAT", zlib.compress(raw, 9)) + chunk(b"IEND", b""))


def overlay(beta0: np.ndarray, layover: np.ndarray, shadow: np.ndarray, valid: np.ndarray,
            lo_db: float = -25.0, hi_db: float = 5.0) -> np.ndarray:
    """beta0 in dB as gray, layover tinted red, shadow tinted blue."""
    db = 10.0 * np.log10(np.maximum(beta0, 1e-6))
    g = np.clip((db - lo_db) / (hi_db - lo_db), 0.0, 1.0) * 255.0
    g = np.where(valid, g, 0.0)
    rgb = np.stack([g, g, g], axis=-1)
    rgb[layover] = 0.6 * rgb[layover] + 0.4 * np.array([255.0, 40.0, 40.0])
    rgb[shadow] = 0.6 * rgb[shadow] + 0.4 * np.array([20.0, 60.0, 200.0])
    return rgb
