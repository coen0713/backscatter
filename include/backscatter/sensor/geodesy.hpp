#pragma once

#include "backscatter/math/vec3.hpp"

namespace bsar {

/// WGS84 geodetic coordinates (degrees, metres above the ellipsoid).
struct Geodetic {
  double lat_deg = 0.0;
  double lon_deg = 0.0;
  double height = 0.0;
};

namespace wgs84 {
inline constexpr double kSemiMajor = 6378137.0;
inline constexpr double kFlattening = 1.0 / 298.257223563;
inline constexpr double kSemiMinor = kSemiMajor * (1.0 - kFlattening);
inline constexpr double kEccentricitySq = kFlattening * (2.0 - kFlattening);
}  // namespace wgs84

Vec3d geodetic_to_ecef(const Geodetic& g);
/// Iterative inverse; converges to sub-millimetre accuracy in a few steps.
Geodetic ecef_to_geodetic(const Vec3d& ecef);

/// Local east-north-up tangent frame anchored at a geodetic origin. The scene
/// lives in this frame; orbits (ECEF) are mapped into it.
class EnuFrame {
 public:
  EnuFrame() = default;
  explicit EnuFrame(const Geodetic& origin);

  [[nodiscard]] const Geodetic& origin() const { return origin_; }
  [[nodiscard]] Vec3d to_enu(const Vec3d& ecef) const;
  [[nodiscard]] Vec3d to_ecef(const Vec3d& enu) const;
  /// Rotate a direction or velocity (no translation).
  [[nodiscard]] Vec3d rotate_to_enu(const Vec3d& v) const;
  [[nodiscard]] Vec3d rotate_to_ecef(const Vec3d& v) const;
  [[nodiscard]] Vec3d geodetic_to_enu(const Geodetic& g) const {
    return to_enu(geodetic_to_ecef(g));
  }

 private:
  Geodetic origin_;
  Vec3d origin_ecef_;
  Vec3d east_{1, 0, 0};
  Vec3d north_{0, 1, 0};
  Vec3d up_{0, 0, 1};
};

}  // namespace bsar
