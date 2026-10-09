#include "backscatter/trace/hillshade.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "backscatter/math/constants.hpp"
#include "backscatter/util/parallel.hpp"

namespace bsar {

Image<float> render_hillshade(const Scene& scene, const HillshadeConfig& config) {
  if (!scene.built()) {
    throw std::invalid_argument("render_hillshade: scene is not built");
  }
  const Aabb& b = scene.bounds();
  const Vec3f ext = b.extent();
  const double ps = config.pixel_size > 0.0 ? config.pixel_size : std::max(ext.x, ext.y) / 512.0;
  const auto w = static_cast<std::size_t>(std::max(1.0, std::floor(ext.x / ps)));
  const auto h = static_cast<std::size_t>(std::max(1.0, std::floor(ext.y / ps)));

  const double az = deg_to_rad(config.sun_azimuth_deg);
  const double el = deg_to_rad(config.sun_elevation_deg);
  const Vec3d sun{std::sin(az) * std::cos(el), std::cos(az) * std::cos(el), std::sin(el)};
  const float top = b.hi.z + 1.0f;

  Image<float> out(w, h);
  parallel_for(0, h, config.threads, [&](std::size_t j) {
    for (std::size_t i = 0; i < w; ++i) {
      const double x = b.lo.x + (static_cast<double>(i) + 0.5) * ps;
      const double y = b.hi.y - (static_cast<double>(j) + 0.5) * ps;
      Ray ray{Vec3f(Vec3d{x, y, top}), Vec3f{0.0f, 0.0f, -1.0f}, 0.0f,
              std::numeric_limits<float>::infinity()};
      Hit hit;
      if (!scene.intersect(ray, hit)) {
        continue;
      }
      Vec3d n(scene.mesh.normal(hit.prim));
      if (n.z < 0.0) {
        n = -n;
      }
      double shade = std::max(0.0, dot(n, sun));
      if (config.cast_shadows && shade > 0.0) {
        const Vec3f p = ray.origin + ray.dir * hit.t;
        const Ray to_sun{p + Vec3f(n) * 1e-3f, Vec3f(sun), 0.0f,
                         std::numeric_limits<float>::infinity()};
        if (scene.occluded(to_sun)) {
          shade = 0.0;
        }
      }
      out.at(i, j) = static_cast<float>(shade);
    }
  });
  return out;
}

}  // namespace bsar
