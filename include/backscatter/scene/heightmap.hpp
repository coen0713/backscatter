#pragma once

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <vector>

#include "backscatter/geometry/mesh.hpp"
#include "backscatter/sensor/geodesy.hpp"

namespace bsar {

/// Regular elevation grid. Sample (i, j) sits at (x0 + i dx, y0 + j dy) with
/// j increasing northwards. For geographic grids x is longitude and y is
/// latitude in degrees; otherwise both are metres in the scene ENU frame.
struct Heightmap {
  std::size_t nx = 0;
  std::size_t ny = 0;
  double x0 = 0.0;
  double y0 = 0.0;
  double dx = 1.0;
  double dy = 1.0;
  bool geographic = false;
  float nodata = -32768.0f;
  std::vector<float> z;  // row-major, j * nx + i

  [[nodiscard]] float at(std::size_t i, std::size_t j) const { return z[j * nx + i]; }
  [[nodiscard]] float& at(std::size_t i, std::size_t j) { return z[j * nx + i]; }
  [[nodiscard]] bool valid(std::size_t i, std::size_t j) const { return at(i, j) != nodata; }
  /// Bilinear height at grid coordinates (x, y); NaN outside or near no-data.
  [[nodiscard]] double sample(double x, double y) const;
};

/// ESRI ASCII grid (.asc), as written by `gdal_translate -of AAIGrid`.
/// Set `geographic` for grids in degrees (e.g. Copernicus GLO-30 tiles).
Heightmap load_esri_ascii(const std::filesystem::path& path, bool geographic);

#ifdef BSAR_HAVE_GDAL
/// Any GDAL-readable raster (GeoTIFF etc.), band 1. Geographic CRS detection
/// uses the dataset's spatial reference.
Heightmap load_dem_gdal(const std::filesystem::path& path);
#endif

/// Two triangles per grid cell, skipping cells that touch no-data. Geographic
/// grids are converted to ENU through `frame`, which is then required.
TriangleMesh heightmap_to_mesh(const Heightmap& hm, std::uint16_t material,
                               const EnuFrame* frame = nullptr);

/// Deterministic synthetic terrain: sum of value-noise octaves.
Heightmap make_fractal_terrain(std::size_t n, double cell_size, double amplitude,
                               std::uint64_t seed);

}  // namespace bsar
