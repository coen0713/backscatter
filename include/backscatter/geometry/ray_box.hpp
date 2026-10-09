#pragma once

#include <algorithm>
#include <cmath>

#include "backscatter/geometry/ray.hpp"
#include "backscatter/math/aabb.hpp"

namespace bsar {

/// Reciprocal ray direction for slab tests. Zero components are replaced by a
/// tiny signed value so that `(plane - origin) * inv` never forms 0 * inf.
inline Vec3f safe_inverse_direction(const Vec3f& d) {
  constexpr float kTiny = 1e-20f;
  auto inv = [](float c) { return 1.0f / (std::abs(c) < kTiny ? std::copysign(kTiny, c) : c); };
  return {inv(d.x), inv(d.y), inv(d.z)};
}

/// Slab test. Returns the entry distance via `t_near` when the box overlaps
/// the ray interval [tmin, tmax].
inline bool intersect_aabb(const Aabb& box, const Vec3f& origin, const Vec3f& inv_dir, float tmin,
                           float tmax, float& t_near) {
  const float tx0 = (box.lo.x - origin.x) * inv_dir.x;
  const float tx1 = (box.hi.x - origin.x) * inv_dir.x;
  const float ty0 = (box.lo.y - origin.y) * inv_dir.y;
  const float ty1 = (box.hi.y - origin.y) * inv_dir.y;
  const float tz0 = (box.lo.z - origin.z) * inv_dir.z;
  const float tz1 = (box.hi.z - origin.z) * inv_dir.z;
  const float t0 = std::max({std::min(tx0, tx1), std::min(ty0, ty1), std::min(tz0, tz1), tmin});
  const float t1 = std::min({std::max(tx0, tx1), std::max(ty0, ty1), std::max(tz0, tz1), tmax});
  t_near = t0;
  return t0 <= t1;
}

}  // namespace bsar
