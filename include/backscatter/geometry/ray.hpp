#pragma once

#include <cstdint>
#include <limits>

#include "backscatter/math/vec3.hpp"

namespace bsar {

inline constexpr std::uint32_t kInvalidPrim = 0xFFFFFFFFu;

/// A ray in scene-local single-precision coordinates. Rays that start at a
/// distant sensor are first clipped to the scene bounds in double precision
/// (see trace/sensor_rays.hpp) so `origin` is always close to the geometry.
struct Ray {
  Vec3f origin;
  Vec3f dir;
  float tmin = 0.0f;
  float tmax = std::numeric_limits<float>::infinity();
};

struct Hit {
  float t = std::numeric_limits<float>::infinity();
  float u = 0.0f;  // barycentric weight of vertex 1
  float v = 0.0f;  // barycentric weight of vertex 2
  std::uint32_t prim = kInvalidPrim;

  [[nodiscard]] bool valid() const { return prim != kInvalidPrim; }
};

}  // namespace bsar
