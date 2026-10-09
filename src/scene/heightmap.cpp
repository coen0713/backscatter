#include "backscatter/scene/heightmap.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>

#include "backscatter/math/rng.hpp"

#ifdef BSAR_HAVE_GDAL
#include <gdal.h>
#include <ogr_srs_api.h>
#endif

namespace bsar {

double Heightmap::sample(double x, double y) const {
  const double fi = (x - x0) / dx;
  const double fj = (y - y0) / dy;
  constexpr double kNan = std::numeric_limits<double>::quiet_NaN();
  if (nx < 2 || ny < 2 || !(fi >= 0.0) || !(fj >= 0.0) || fi > static_cast<double>(nx - 1) ||
      fj > static_cast<double>(ny - 1)) {
    return kNan;
  }
  const auto i = std::min(static_cast<std::size_t>(fi), nx - 2);
  const auto j = std::min(static_cast<std::size_t>(fj), ny - 2);
  if (!valid(i, j) || !valid(i + 1, j) || !valid(i, j + 1) || !valid(i + 1, j + 1)) {
    return kNan;
  }
  const double u = fi - static_cast<double>(i);
  const double v = fj - static_cast<double>(j);
  return (1 - u) * (1 - v) * at(i, j) + u * (1 - v) * at(i + 1, j) + (1 - u) * v * at(i, j + 1) +
         u * v * at(i + 1, j + 1);
}

Heightmap load_esri_ascii(const std::filesystem::path& path, bool geographic) {
  std::ifstream in(path);
  if (!in) {
    throw std::runtime_error("load_esri_ascii: cannot open " + path.string());
  }
  Heightmap hm;
  hm.geographic = geographic;
  double xll = 0.0;
  double yll = 0.0;
  double cell_x = 0.0;
  double cell_y = 0.0;
  bool corner = true;
  std::size_t rows = 0;
  std::size_t cols = 0;

  // Header: "key value" lines until the first numeric token.
  std::string key;
  while (in >> key) {
    std::string lower = key;
    std::transform(lower.begin(), lower.end(), lower.begin(),
                   [](unsigned char ch) { return static_cast<char>(std::tolower(ch)); });
    if (!lower.empty() && (std::isdigit(static_cast<unsigned char>(lower[0])) != 0 ||
                           lower[0] == '-' || lower[0] == '.')) {
      break;  // first data value consumed into `key`
    }
    double value = 0.0;
    in >> value;
    if (lower == "ncols") {
      cols = static_cast<std::size_t>(value);
    } else if (lower == "nrows") {
      rows = static_cast<std::size_t>(value);
    } else if (lower == "xllcorner" || lower == "xllcenter") {
      xll = value;
      corner = lower == "xllcorner";
    } else if (lower == "yllcorner" || lower == "yllcenter") {
      yll = value;
    } else if (lower == "cellsize") {
      cell_x = cell_y = value;
    } else if (lower == "dx") {
      cell_x = value;
    } else if (lower == "dy") {
      cell_y = value;
    } else if (lower == "nodata_value") {
      hm.nodata = static_cast<float>(value);
    } else {
      throw std::runtime_error("load_esri_ascii: unknown header key '" + key + "'");
    }
  }
  if (rows < 2 || cols < 2 || cell_x <= 0.0 || cell_y <= 0.0) {
    throw std::runtime_error("load_esri_ascii: incomplete header in " + path.string());
  }
  hm.nx = cols;
  hm.ny = rows;
  hm.dx = cell_x;
  hm.dy = cell_y;
  hm.x0 = corner ? xll + 0.5 * cell_x : xll;
  hm.y0 = corner ? yll + 0.5 * cell_y : yll;
  hm.z.resize(rows * cols);

  std::size_t k = 0;
  auto store = [&](double value) {
    const std::size_t r = k / cols;  // row from the top (north)
    const std::size_t c = k % cols;
    hm.at(c, rows - 1 - r) = static_cast<float>(value);
    ++k;
  };
  store(std::stod(key));
  double value = 0.0;
  while (k < rows * cols && in >> value) {
    store(value);
  }
  if (k != rows * cols) {
    throw std::runtime_error("load_esri_ascii: expected " + std::to_string(rows * cols) +
                             " values, read " + std::to_string(k));
  }
  return hm;
}

#ifdef BSAR_HAVE_GDAL
Heightmap load_dem_gdal(const std::filesystem::path& path) {
  GDALAllRegister();
  GDALDatasetH ds = GDALOpen(path.string().c_str(), GA_ReadOnly);
  if (ds == nullptr) {
    throw std::runtime_error("load_dem_gdal: cannot open " + path.string());
  }
  double gt[6];
  if (GDALGetGeoTransform(ds, gt) != CE_None || gt[2] != 0.0 || gt[4] != 0.0) {
    GDALClose(ds);
    throw std::runtime_error("load_dem_gdal: need a north-up geotransform");
  }
  Heightmap hm;
  hm.nx = static_cast<std::size_t>(GDALGetRasterXSize(ds));
  hm.ny = static_cast<std::size_t>(GDALGetRasterYSize(ds));
  hm.dx = gt[1];
  hm.dy = -gt[5];
  hm.x0 = gt[0] + 0.5 * gt[1];
  hm.y0 = gt[3] + (static_cast<double>(hm.ny) - 0.5) * gt[5];

  OGRSpatialReferenceH srs = OSRNewSpatialReference(GDALGetProjectionRef(ds));
  hm.geographic = srs != nullptr && OSRIsGeographic(srs) != 0;
  if (srs != nullptr) {
    OSRDestroySpatialReference(srs);
  }

  GDALRasterBandH band = GDALGetRasterBand(ds, 1);
  int has_nodata = 0;
  const double nodata = GDALGetRasterNoDataValue(band, &has_nodata);
  if (has_nodata != 0) {
    hm.nodata = static_cast<float>(nodata);
  }
  std::vector<float> rows(hm.nx * hm.ny);
  if (GDALRasterIO(band, GF_Read, 0, 0, static_cast<int>(hm.nx), static_cast<int>(hm.ny),
                   rows.data(), static_cast<int>(hm.nx), static_cast<int>(hm.ny), GDT_Float32, 0,
                   0) != CE_None) {
    GDALClose(ds);
    throw std::runtime_error("load_dem_gdal: read failed for " + path.string());
  }
  GDALClose(ds);
  hm.z.resize(rows.size());
  for (std::size_t r = 0; r < hm.ny; ++r) {
    std::copy_n(rows.begin() + static_cast<std::ptrdiff_t>(r * hm.nx), hm.nx,
                hm.z.begin() + static_cast<std::ptrdiff_t>((hm.ny - 1 - r) * hm.nx));
  }
  return hm;
}
#endif

TriangleMesh heightmap_to_mesh(const Heightmap& hm, std::uint16_t material, const EnuFrame* frame) {
  if (hm.geographic && frame == nullptr) {
    throw std::invalid_argument("heightmap_to_mesh: geographic grid needs an ENU frame");
  }
  TriangleMesh mesh;
  if (hm.nx < 2 || hm.ny < 2) {
    return mesh;
  }
  mesh.reserve(hm.nx * hm.ny, 2 * (hm.nx - 1) * (hm.ny - 1));
  for (std::size_t j = 0; j < hm.ny; ++j) {
    for (std::size_t i = 0; i < hm.nx; ++i) {
      const double x = hm.x0 + static_cast<double>(i) * hm.dx;
      const double y = hm.y0 + static_cast<double>(j) * hm.dy;
      const double z = hm.valid(i, j) ? hm.at(i, j) : 0.0;
      const Vec3d p = hm.geographic ? frame->geodetic_to_enu({y, x, z}) : Vec3d{x, y, z};
      mesh.add_vertex(Vec3f(p));
    }
  }
  for (std::size_t j = 0; j + 1 < hm.ny; ++j) {
    for (std::size_t i = 0; i + 1 < hm.nx; ++i) {
      if (!hm.valid(i, j) || !hm.valid(i + 1, j) || !hm.valid(i, j + 1) ||
          !hm.valid(i + 1, j + 1)) {
        continue;
      }
      const auto v00 = static_cast<std::uint32_t>(j * hm.nx + i);
      const auto v10 = v00 + 1;
      const auto v01 = static_cast<std::uint32_t>(v00 + hm.nx);
      const auto v11 = v01 + 1;
      mesh.add_triangle(v00, v10, v11, material);
      mesh.add_triangle(v00, v11, v01, material);
    }
  }
  return mesh;
}

Heightmap make_fractal_terrain(std::size_t n, double cell_size, double amplitude,
                               std::uint64_t seed) {
  Heightmap hm;
  hm.nx = hm.ny = n;
  hm.dx = hm.dy = cell_size;
  hm.x0 = hm.y0 = -0.5 * static_cast<double>(n - 1) * cell_size;
  hm.z.assign(n * n, 0.0f);
  const CounterRng rng(seed);

  auto lattice = [&](std::uint64_t octave, std::int64_t ix, std::int64_t iy) {
    const auto key =
        (static_cast<std::uint64_t>(ix) & 0xFFFFFFFFull) | (static_cast<std::uint64_t>(iy) << 32);
    return rng.uniform4(octave, key)[0] * 2.0 - 1.0;
  };
  auto smooth = [](double t) { return t * t * (3.0 - 2.0 * t); };

  constexpr int kOctaves = 7;
  double freq = 3.0 / static_cast<double>(n);  // ~3 features across the tile
  double amp = 1.0;
  double total = 0.0;
  for (int o = 0; o < kOctaves; ++o) {
    for (std::size_t j = 0; j < n; ++j) {
      for (std::size_t i = 0; i < n; ++i) {
        const double fx = static_cast<double>(i) * freq;
        const double fy = static_cast<double>(j) * freq;
        const auto ix = static_cast<std::int64_t>(std::floor(fx));
        const auto iy = static_cast<std::int64_t>(std::floor(fy));
        const double u = smooth(fx - static_cast<double>(ix));
        const double v = smooth(fy - static_cast<double>(iy));
        const auto oo = static_cast<std::uint64_t>(o);
        const double a = lattice(oo, ix, iy);
        const double b = lattice(oo, ix + 1, iy);
        const double c = lattice(oo, ix, iy + 1);
        const double d = lattice(oo, ix + 1, iy + 1);
        const double value = (1 - u) * (1 - v) * a + u * (1 - v) * b + (1 - u) * v * c + u * v * d;
        hm.at(i, j) += static_cast<float>(amp * value);
      }
    }
    total += amp;
    freq *= 2.0;
    amp *= 0.5;
  }
  const auto scale = static_cast<float>(amplitude / total);
  for (auto& z : hm.z) {
    z *= scale;
  }
  return hm;
}

}  // namespace bsar
