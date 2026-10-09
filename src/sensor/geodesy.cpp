#include "backscatter/sensor/geodesy.hpp"

#include <cmath>

#include "backscatter/math/constants.hpp"

namespace bsar {

Vec3d geodetic_to_ecef(const Geodetic& g) {
  const double lat = deg_to_rad(g.lat_deg);
  const double lon = deg_to_rad(g.lon_deg);
  const double s = std::sin(lat);
  const double n = wgs84::kSemiMajor / std::sqrt(1.0 - wgs84::kEccentricitySq * s * s);
  const double r = (n + g.height) * std::cos(lat);
  return {r * std::cos(lon), r * std::sin(lon),
          (n * (1.0 - wgs84::kEccentricitySq) + g.height) * s};
}

Geodetic ecef_to_geodetic(const Vec3d& p) {
  const double lon = std::atan2(p.y, p.x);
  const double rho = std::hypot(p.x, p.y);
  double lat = std::atan2(p.z, rho * (1.0 - wgs84::kEccentricitySq));
  double h = 0.0;
  for (int i = 0; i < 8; ++i) {
    const double s = std::sin(lat);
    const double n = wgs84::kSemiMajor / std::sqrt(1.0 - wgs84::kEccentricitySq * s * s);
    // Use whichever form is well conditioned at this latitude.
    if (std::abs(std::cos(lat)) > 1e-3) {
      h = rho / std::cos(lat) - n;
    } else {
      h = p.z / s - n * (1.0 - wgs84::kEccentricitySq);
    }
    lat = std::atan2(p.z, rho * (1.0 - wgs84::kEccentricitySq * n / (n + h)));
  }
  return {rad_to_deg(lat), rad_to_deg(lon), h};
}

EnuFrame::EnuFrame(const Geodetic& origin)
    : origin_(origin), origin_ecef_(geodetic_to_ecef(origin)) {
  const double lat = deg_to_rad(origin.lat_deg);
  const double lon = deg_to_rad(origin.lon_deg);
  east_ = {-std::sin(lon), std::cos(lon), 0.0};
  north_ = {-std::sin(lat) * std::cos(lon), -std::sin(lat) * std::sin(lon), std::cos(lat)};
  up_ = {std::cos(lat) * std::cos(lon), std::cos(lat) * std::sin(lon), std::sin(lat)};
}

Vec3d EnuFrame::rotate_to_enu(const Vec3d& v) const {
  return {dot(east_, v), dot(north_, v), dot(up_, v)};
}

Vec3d EnuFrame::rotate_to_ecef(const Vec3d& v) const {
  return east_ * v.x + north_ * v.y + up_ * v.z;
}

Vec3d EnuFrame::to_enu(const Vec3d& ecef) const { return rotate_to_enu(ecef - origin_ecef_); }

Vec3d EnuFrame::to_ecef(const Vec3d& enu) const { return origin_ecef_ + rotate_to_ecef(enu); }

}  // namespace bsar
