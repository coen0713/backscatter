#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>

#include "backscatter/geometry/ray.hpp"
#include "backscatter/math/aabb.hpp"
#include "backscatter/math/vec3.hpp"

namespace bsar {

struct ClippedRay {
  Ray ray;          // scene-local float ray starting just outside the bounds
  double t_offset;  // distance from the true origin to ray.origin
};

/// Clip a ray from a possibly distant origin (a satellite ~10^6 m away) to
/// the scene bounds in double precision, then hand back a float ray that
/// starts close to the geometry. Total path length = t_offset + hit.t.
inline std::optional<ClippedRay> clip_to_bounds(const Vec3d& origin, const Vec3d& dir,
                                                const Aabb& bounds, double margin) {
  if (bounds.empty()) {
    return std::nullopt;
  }
  const Vec3d lo = Vec3d(bounds.lo) - Vec3d{margin, margin, margin};
  const Vec3d hi = Vec3d(bounds.hi) + Vec3d{margin, margin, margin};
  double t0 = 0.0;
  double t1 = std::numeric_limits<double>::infinity();
  for (std::size_t a = 0; a < 3; ++a) {
    if (std::abs(dir[a]) < 1e-300) {
      if (origin[a] < lo[a] || origin[a] > hi[a]) {
        return std::nullopt;
      }
      continue;
    }
    const double inv = 1.0 / dir[a];
    double ta = (lo[a] - origin[a]) * inv;
    double tb = (hi[a] - origin[a]) * inv;
    if (ta > tb) {
      std::swap(ta, tb);
    }
    t0 = std::max(t0, ta);
    t1 = std::min(t1, tb);
  }
  if (t0 > t1) {
    return std::nullopt;
  }
  ClippedRay out;
  out.t_offset = t0;
  out.ray.origin = Vec3f(origin + dir * t0);
  out.ray.dir = Vec3f(dir);
  out.ray.tmin = 0.0f;
  out.ray.tmax = static_cast<float>(t1 - t0) + 1.0f;
  return out;
}

}  // namespace bsar
