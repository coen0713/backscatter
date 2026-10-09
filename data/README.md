# Data

No data is checked in. These scripts fetch everything from public sources into
`data/raw/` and `data/cache/` (both git-ignored).

| Script | Source | Account needed |
|---|---|---|
| `fetch_dem.py` | Copernicus GLO-30 DEM, AWS Open Data | no |
| `fetch_osm_buildings.py` | OpenStreetMap via the Overpass API | no |
| `fetch_orbit.py` | Sentinel-1 POEORB precise orbits, ESA SNAP mirror | no |

Sentinel-1 GRD scenes for the real-data comparison (`docs/validation.md`,
section 2) come from the [Copernicus Data Space Ecosystem](https://dataspace.copernicus.eu/),
which needs a free account. Download them by hand; nothing here automates logging in.

DEM conversion uses the GDAL command-line tools. Without GDAL, any ESRI ASCII
grid (`.asc`) works directly.

## End-to-end example

```bash
python data/fetch_dem.py --lat 46 --lon 7 --bbox 7.95 46.53 8.05 46.60 --ellipsoidal --out data/raw/eiger.asc
python data/fetch_orbit.py --mission S1A --time 2024-06-01T05:30:00
build/release/apps/bsar_render --scene dem:data/raw/eiger.asc --geographic \
    --orbit data/raw/S1A_OPER_AUX_POEORB_<...>.EOF --time 2024-06-01T05:30:00 \
    --range-spacing 10 --azimuth-spacing 10 --out renders/eiger
```

Footprint file format (written by `fetch_osm_buildings.py`):

```
crs wgs84            # or "crs enu" for scene-local metres
building 24.0 concrete
7.95112 46.53321     # lon lat (or x y)
...
end
```
