#pragma once

#include <cmath>
#include <utility>

#include "backscatter/geometry/ray.hpp"
#include "backscatter/math/vec3.hpp"

namespace bsar {

/// Per-ray precomputation for the watertight ray/triangle test of Woop, Benthin
/// and Wald, "Watertight Ray/Triangle Intersection" (JCGT 2013).
///
/// The ray is sheared and scaled so it points along +z; the test then reduces
/// to 2D edge functions that are evaluated consistently for edges shared by
/// neighbouring triangles. A ray can therefore never slip between two
/// triangles that share an edge, which matters for SAR: terrain is viewed at
/// grazing angles where conventional tests leak through mesh seams.
///
/// The guarantee needs `a*b - c*d` to round both products separately, so code
/// including this header must be compiled without FMA contraction
/// (-ffp-contract=off, set as a usage requirement of the library target).
struct WatertightRay {
  Vec3f org;
  int kx = 0;
  int ky = 1;
  int kz = 2;
  float sx = 0.0f;
  float sy = 0.0f;
  float sz = 1.0f;

  explicit WatertightRay(const Ray& ray) : org(ray.origin) {
    kz = max_dimension(abs(ray.dir));
    kx = (kz + 1) % 3;
    ky = (kx + 1) % 3;
    if (ray.dir[static_cast<std::size_t>(kz)] < 0.0f) {
      std::swap(kx, ky);
    }
    const float dz = ray.dir[static_cast<std::size_t>(kz)];
    sx = ray.dir[static_cast<std::size_t>(kx)] / dz;
    sy = ray.dir[static_cast<std::size_t>(ky)] / dz;
    sz = 1.0f / dz;
  }
};

/// Returns true for a hit with tmin < t < tmax and writes t and the
/// barycentric weights (u for v1, v for v2). Back faces are reported too.
inline bool intersect_triangle(const WatertightRay& r, const Vec3f& v0, const Vec3f& v1,
                               const Vec3f& v2, float tmin, float tmax, float& t_out, float& u_out,
                               float& v_out) {
  const auto kx = static_cast<std::size_t>(r.kx);
  const auto ky = static_cast<std::size_t>(r.ky);
  const auto kz = static_cast<std::size_t>(r.kz);

  const Vec3f a = v0 - r.org;
  const Vec3f b = v1 - r.org;
  const Vec3f c = v2 - r.org;

  const float ax = a[kx] - r.sx * a[kz];
  const float ay = a[ky] - r.sy * a[kz];
  const float bx = b[kx] - r.sx * b[kz];
  const float by = b[ky] - r.sy * b[kz];
  const float cx = c[kx] - r.sx * c[kz];
  const float cy = c[ky] - r.sy * c[kz];

  float u = cx * by - cy * bx;
  float v = ax * cy - ay * cx;
  float w = bx * ay - by * ax;

  // Fall back to double precision when an edge function is exactly zero, so
  // the sign of an edge test is never decided by rounding.
  if (u == 0.0f || v == 0.0f || w == 0.0f) {
    const double cxby = static_cast<double>(cx) * static_cast<double>(by);
    const double cybx = static_cast<double>(cy) * static_cast<double>(bx);
    u = static_cast<float>(cxby - cybx);
    const double axcy = static_cast<double>(ax) * static_cast<double>(cy);
    const double aycx = static_cast<double>(ay) * static_cast<double>(cx);
    v = static_cast<float>(axcy - aycx);
    const double bxay = static_cast<double>(bx) * static_cast<double>(ay);
    const double byax = static_cast<double>(by) * static_cast<double>(ax);
    w = static_cast<float>(bxay - byax);
  }

  if ((u < 0.0f || v < 0.0f || w < 0.0f) && (u > 0.0f || v > 0.0f || w > 0.0f)) {
    return false;
  }
  const float det = u + v + w;
  if (det == 0.0f) {
    return false;
  }

  const float az = r.sz * a[kz];
  const float bz = r.sz * b[kz];
  const float cz = r.sz * c[kz];
  const float t_scaled = u * az + v * bz + w * cz;

  const float inv_det = 1.0f / det;
  const float t = t_scaled * inv_det;
  if (!(t > tmin && t < tmax)) {
    return false;
  }
  t_out = t;
  u_out = v * inv_det;
  v_out = w * inv_det;
  return true;
}

}  // namespace bsar
