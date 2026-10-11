"""Inputs on common grids: ellipsoidal DEM crops, WorldCover, UTM tile grids.

Uses rasterio (pip wheels bundle GDAL, so this also works on Windows). Source
files are downloaded whole into data/cache once (see fetch.py) and read
locally.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import rasterio
from rasterio.crs import CRS
from rasterio.enums import WktVersion
from rasterio.enums import Resampling
from rasterio.merge import merge
from rasterio.transform import Affine
from rasterio.warp import reproject, transform_bounds
from rasterio.windows import from_bounds

from validation.fetch import download

DEM_URL = ("https://copernicus-dem-30m.s3.amazonaws.com/Copernicus_DSM_COG_10_{ns}{lat:02d}_00_"
           "{ew}{lon:03d}_00_DEM/Copernicus_DSM_COG_10_{ns}{lat:02d}_00_{ew}{lon:03d}_00_DEM.tif")
GEOID_URL = "https://cdn.proj.org/us_nga_egm08_25.tif"  # EGM2008 undulation N, 2.5 arcmin
WORLDCOVER_URL = ("https://esa-worldcover.s3.eu-central-1.amazonaws.com/v200/2021/map/"
                  "ESA_WorldCover_10m_2021_v200_{ns}{lat:02d}{ew}{lon:03d}_Map.tif")

WGS84 = CRS.from_epsg(4326)
UTM32N = CRS.from_epsg(32632)
NODATA = -9999.0


@dataclass
class Grid:
    """A north-up raster grid."""

    transform: Affine
    width: int
    height: int
    crs: CRS

    @property
    def bounds(self) -> tuple[float, float, float, float]:
        t = self.transform
        return (t.c, t.f + t.e * self.height, t.c + t.a * self.width, t.f)


def _dem_tile(lat: int, lon: int) -> Path:
    url = DEM_URL.format(ns="N" if lat >= 0 else "S", lat=abs(lat), ew="E" if lon >= 0 else "W",
                         lon=abs(lon))
    return download(url)


def ellipsoidal_dem(bounds: tuple[float, float, float, float]) -> tuple[np.ndarray, Grid]:
    """Copernicus GLO-30 heights over `bounds` (lon/lat), converted from
    EGM2008 geoid heights to WGS84 ellipsoidal heights (h = H + N), which is
    what orbit geometry needs. Returns (heights, grid); row 0 is north."""
    w, s, e, n = bounds
    paths = [_dem_tile(lat, lon)
             for lat in range(math.floor(s), math.floor(n - 1e-9) + 1)
             for lon in range(math.floor(w), math.floor(e - 1e-9) + 1)]
    sources = [rasterio.open(p) for p in paths]
    try:
        mosaic, transform = merge(sources, bounds=bounds, nodata=NODATA)
    finally:
        for src in sources:
            src.close()
    heights = mosaic[0].astype(np.float64)
    grid = Grid(transform, heights.shape[1], heights.shape[0], WGS84)

    undulation = np.zeros_like(heights)
    with rasterio.open(download(GEOID_URL)) as geoid:
        reproject(rasterio.band(geoid, 1), undulation, dst_transform=transform, dst_crs=WGS84,
                  resampling=Resampling.bilinear)
    valid = heights != NODATA
    heights[valid] += undulation[valid]
    return heights.astype(np.float32), grid


def write_asc(path: Path, values: np.ndarray, grid: Grid, nodata: float = NODATA) -> None:
    """ESRI ASCII grid + .prj, readable by bsar_render and by GDAL."""
    t = grid.transform
    if t.b != 0 or t.d != 0:
        raise ValueError("write_asc: grid must be north-up")
    with open(path, "w") as f:
        f.write(f"ncols {grid.width}\nnrows {grid.height}\n")
        f.write(f"xllcorner {t.c:.12f}\nyllcorner {t.f + t.e * grid.height:.12f}\n")
        if abs(t.a) == abs(t.e):
            f.write(f"cellsize {t.a:.15g}\n")
        else:
            f.write(f"dx {t.a:.15g}\ndy {-t.e:.15g}\n")
        f.write(f"NODATA_value {nodata:g}\n")
        np.savetxt(f, np.where(np.isnan(values), nodata, values), fmt="%.3f")
    path.with_suffix(".prj").write_text(grid.crs.to_wkt(version=WktVersion.WKT1_ESRI))


def read_raster(path: Path) -> tuple[np.ndarray, Grid, float | None]:
    with rasterio.open(path) as ds:
        crs = ds.crs or WGS84
        return ds.read(1), Grid(ds.transform, ds.width, ds.height, crs), ds.nodata


def utm_grid(bounds_lonlat: tuple[float, float, float, float], resolution: float = 30.0,
             crs: CRS = UTM32N) -> Grid:
    """The tile's footprint in UTM, snapped outwards to multiples of `resolution`."""
    x0, y0, x1, y1 = transform_bounds(WGS84, crs, *bounds_lonlat, densify_pts=21)
    x0 = math.floor(x0 / resolution) * resolution
    y0 = math.floor(y0 / resolution) * resolution
    x1 = math.ceil(x1 / resolution) * resolution
    y1 = math.ceil(y1 / resolution) * resolution
    return Grid(Affine(resolution, 0, x0, 0, -resolution, y1), round((x1 - x0) / resolution),
                round((y1 - y0) / resolution), crs)


def resample(values: np.ndarray, src: Grid, dst: Grid, resampling: Resampling,
             src_nodata: float | None, dst_nodata: float = NODATA) -> np.ndarray:
    out = np.full((dst.height, dst.width), dst_nodata, dtype=np.float64)
    reproject(values.astype(np.float64), out, src_transform=src.transform, src_crs=src.crs,
              src_nodata=src_nodata, dst_transform=dst.transform, dst_crs=dst.crs,
              dst_nodata=dst_nodata, resampling=resampling)
    return out


def worldcover(dst: Grid) -> np.ndarray:
    """ESA WorldCover 2021 class codes on `dst`, by majority (mode) of the 10 m
    pixels. 0 = no data."""
    w, s, e, n = transform_bounds(dst.crs, WGS84, *dst.bounds, densify_pts=21)
    pad = 0.01
    out = np.zeros((dst.height, dst.width), dtype=np.uint8)
    # WorldCover comes in 3 x 3 degree tiles; paint each overlapping one in turn.
    for lat3 in range(math.floor((s - pad) / 3) * 3, math.floor((n + pad) / 3) * 3 + 1, 3):
        for lon3 in range(math.floor((w - pad) / 3) * 3, math.floor((e + pad) / 3) * 3 + 1, 3):
            url = WORLDCOVER_URL.format(ns="N" if lat3 >= 0 else "S", lat=abs(lat3),
                                        ew="E" if lon3 >= 0 else "W", lon=abs(lon3))
            with rasterio.open(download(url)) as ds:
                bw, bs = max(w - pad, lon3), max(s - pad, lat3)
                be, bn = min(e + pad, lon3 + 3), min(n + pad, lat3 + 3)
                window = from_bounds(bw, bs, be, bn, ds.transform)
                data = ds.read(1, window=window)
                reproject(data, out, src_transform=ds.window_transform(window), src_crs=ds.crs,
                          src_nodata=0, dst_transform=dst.transform, dst_crs=dst.crs,
                          dst_nodata=0, init_dest_nodata=False, resampling=Resampling.mode)
    return out
