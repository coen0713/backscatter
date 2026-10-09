#pragma once

#include <limits>

#include "backscatter/math/vec3.hpp"

namespace bsar {

/// Axis-aligned bounding box. Default-constructed boxes are empty (lo > hi).
struct Aabb {
  Vec3f lo{std::numeric_limits<float>::infinity(), std::numeric_limits<float>::infinity(),
           std::numeric_limits<float>::infinity()};
  Vec3f hi{-std::numeric_limits<float>::infinity(), -std::numeric_limits<float>::infinity(),
           -std::numeric_limits<float>::infinity()};

  constexpr void expand(const Vec3f& p) {
    lo = min(lo, p);
    hi = max(hi, p);
  }
  constexpr void expand(const Aabb& b) {
    lo = min(lo, b.lo);
    hi = max(hi, b.hi);
  }

  [[nodiscard]] constexpr bool empty() const { return lo.x > hi.x || lo.y > hi.y || lo.z > hi.z; }
  [[nodiscard]] constexpr Vec3f extent() const { return empty() ? Vec3f{} : hi - lo; }
  [[nodiscard]] constexpr Vec3f center() const { return (lo + hi) * 0.5f; }

  [[nodiscard]] constexpr float surface_area() const {
    if (empty()) {
      return 0.0f;
    }
    const Vec3f e = hi - lo;
    return 2.0f * (e.x * e.y + e.y * e.z + e.z * e.x);
  }

  [[nodiscard]] constexpr int largest_axis() const { return max_dimension(extent()); }

  /// Corner `i` in [0, 8): bit 0 selects x, bit 1 y, bit 2 z.
  [[nodiscard]] constexpr Vec3f corner(int i) const {
    return {(i & 1) != 0 ? hi.x : lo.x, (i & 2) != 0 ? hi.y : lo.y, (i & 4) != 0 ? hi.z : lo.z};
  }
};

}  // namespace bsar
