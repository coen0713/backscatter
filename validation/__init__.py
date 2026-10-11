"""Real-data validation pipeline (docs/prereg.md).

Modules:
  fetch    - resumable, retrying downloads into data/cache
  regions  - tiles, blocks and the calibration / held-out split
  prepare  - DEM (ellipsoidal heights), WorldCover, common UTM grids
  metrics  - the pre-registered metrics, pure NumPy
  pipeline - per-tile simulation and comparison
  dry_run  - end-to-end run on calibration tiles with self-consistency checks
"""
