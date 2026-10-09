#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#include "backscatter/geometry/mesh.hpp"
#include "backscatter/scene/heightmap.hpp"
#include "backscatter/sensor/geodesy.hpp"

namespace bsar {

struct Vec2d {
  double x = 0.0;
  double y = 0.0;
};

/// A building footprint to extrude. Coordinates are scene ENU metres, or
/// (lon, lat) degrees when loaded from a WGS84 file.
struct Footprint {
  std::vector<Vec2d> ring;  // closed implicitly; either winding is accepted
  double height = 10.0;     // metres above the base
  std::string material = "concrete";
};

struct FootprintSet {
  bool geographic = false;
  std::vector<Footprint> footprints;
};

/// Text format written by data/fetch_osm_buildings.py:
///
///   # comment
///   crs wgs84            (or "crs enu")
///   building <height_m> [material]
///   <x> <y>
///   ...
///   end
FootprintSet load_footprints(const std::filesystem::path& path);

/// Extrude footprints into closed prisms (walls + flat roof; the floor is
/// omitted since it is never visible). The base sits at the lowest terrain
/// height under the footprint so walls reach into the ground.
/// `material_index` maps a material name to the scene's material id.
template <typename MaterialIndex>
void extrude_footprints(const FootprintSet& set, const Heightmap* terrain, const EnuFrame* frame,
                        MaterialIndex&& material_index, TriangleMesh& out);

/// Non-template core used by extrude_footprints.
void extrude_footprint(std::span<const Vec2d> ring_enu, double base, double top,
                       std::uint16_t material, TriangleMesh& out);

/// Ear-clipping triangulation of a simple polygon. Returns index triples
/// (counter-clockwise).
std::vector<std::uint32_t> triangulate_polygon(std::span<const Vec2d> ring);

/// Signed area (positive for counter-clockwise rings).
double signed_area(std::span<const Vec2d> ring);

template <typename MaterialIndex>
void extrude_footprints(const FootprintSet& set, const Heightmap* terrain, const EnuFrame* frame,
                        MaterialIndex&& material_index, TriangleMesh& out) {
  if (set.geographic && frame == nullptr) {
    throw std::invalid_argument("extrude_footprints: geographic footprints need an ENU frame");
  }
  for (const Footprint& fp : set.footprints) {
    if (fp.ring.size() < 3) {
      continue;
    }
    // Terrain grids and footprints share a CRS (both geographic or both ENU).
    double base = std::numeric_limits<double>::infinity();
    if (terrain != nullptr) {
      for (const Vec2d& p : fp.ring) {
        const double h = terrain->sample(p.x, p.y);
        if (!std::isnan(h)) {
          base = std::min(base, h);
        }
      }
    }
    if (!(base < std::numeric_limits<double>::infinity())) {
      base = 0.0;
    }
    std::vector<Vec2d> enu;
    enu.reserve(fp.ring.size());
    double base_z = set.geographic ? std::numeric_limits<double>::infinity() : base;
    for (const Vec2d& p : fp.ring) {
      if (set.geographic) {
        // Earth curvature lowers ENU z away from the origin; take the lowest.
        const Vec3d e = frame->geodetic_to_enu({p.y, p.x, base});
        enu.push_back({e.x, e.y});
        base_z = std::min(base_z, e.z);
      } else {
        enu.push_back(p);
      }
    }
    extrude_footprint(enu, base_z, base_z + fp.height, material_index(fp.material), out);
  }
}

}  // namespace bsar
