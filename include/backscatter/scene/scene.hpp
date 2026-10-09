#pragma once

#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include "backscatter/geometry/bvh.hpp"
#include "backscatter/geometry/mesh.hpp"
#include "backscatter/geometry/wide_bvh.hpp"
#include "backscatter/scene/material.hpp"
#include "backscatter/sensor/geodesy.hpp"

namespace bsar {

enum class AccelKind { Binary, Wide4 };

/// Triangles, materials, an optional geodetic anchor and the acceleration
/// structure. Call build() after editing `mesh`.
class Scene {
 public:
  TriangleMesh mesh;
  std::vector<Material> materials;
  std::optional<EnuFrame> frame;  // set for scenes anchored on the Earth

  /// Id of the material with this name, adding the preset if necessary.
  std::uint16_t material_id(const std::string& name);
  std::uint16_t add_material(const Material& m);

  void build(const BvhBuildConfig& config = {}, AccelKind accel = AccelKind::Wide4);
  [[nodiscard]] bool built() const { return built_; }

  bool intersect(const Ray& ray, Hit& hit) const;
  [[nodiscard]] bool occluded(const Ray& ray) const;

  [[nodiscard]] const Bvh& bvh() const { return bvh_; }
  [[nodiscard]] const Bvh4& bvh4() const { return bvh4_; }
  [[nodiscard]] const Aabb& bounds() const { return bounds_; }
  [[nodiscard]] const Material& material_of(std::uint32_t prim) const {
    return materials[mesh.material(prim)];
  }

 private:
  Bvh bvh_;
  Bvh4 bvh4_;
  AccelKind accel_ = AccelKind::Wide4;
  Aabb bounds_;
  bool built_ = false;
};

/// Flat square of side `size` centred on the origin at z = 0.
Scene make_plane_scene(double size, const std::string& material = "dry_soil");

/// A ridge running north-south through the origin on flat ground. The west
/// face rises at `west_slope_deg`, the east face falls at `east_slope_deg`.
/// Good for checking layover and shadow against closed-form angles.
Scene make_ridge_scene(double west_slope_deg, double east_slope_deg, double height,
                       double ground_size = 2000.0);

/// Flat ground with a single box building centred at the origin: the
/// canonical dihedral (wall-ground) double-bounce target.
Scene make_dihedral_scene(double ground_size, double width, double depth, double height);

/// Named synthetic scenes for the CLI: plane, ridge, dihedral, city, mountains.
Scene make_synthetic_scene(std::string_view name);

}  // namespace bsar
