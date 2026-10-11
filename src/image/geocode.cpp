#include "backscatter/image/geocode.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "backscatter/sensor/sar_geometry.hpp"
#include "backscatter/util/parallel.hpp"

namespace bsar {

GeocodedProducts geocode(const GeometricImage& image, const Trajectory& trajectory,
                         const Scene& scene, std::span<const Vec3d> points, std::size_t width,
                         std::size_t height, unsigned threads) {
  if (points.size() != width * height) {
    throw std::invalid_argument("geocode: points must have width * height entries");
  }
  GeocodedProducts out;
  out.beta0 = Image<float>(width, height);
  out.layover = Image<std::uint8_t>(width, height);
  out.shadow = Image<std::uint8_t>(width, height);
  out.valid = Image<std::uint8_t>(width, height);
  out.range = Image<float>(width, height);
  const std::size_t lines = image.intensity.height;
  const std::size_t bins = image.intensity.width;
  if (lines < 2 || bins < 2) {
    return out;
  }

  parallel_for(0, height, threads, [&](std::size_t j) {
    bool have_guess = false;
    double guess = 0.0;
    for (std::size_t i = 0; i < width; ++i) {
      const Vec3d& x = points[j * width + i];
      if (std::isnan(x.z)) {
        continue;
      }
      // Warm start from the neighbouring cell; a coarse search for the first.
      const double t0 = have_guess ? guess : find_zero_doppler_time(trajectory, x);
      const ZeroDopplerSolution zd = solve_zero_doppler(trajectory, x, t0);
      if (!zd.converged) {
        have_guess = false;
        continue;
      }
      guess = zd.t;
      have_guess = true;
      const double fl = (zd.t - image.t_start) / image.line_interval;
      const double fb = (zd.range - image.near_range) / image.range_spacing - 0.5;
      out.range.at(i, j) = static_cast<float>(zd.range);
      // Ground shadow: start 1 m along the line of sight so the surface the
      // point sits on is not hit at t ~ 0.
      const Vec3d to_sensor = normalize(trajectory.position(zd.t) - x);
      const Ray ray{Vec3f(x), Vec3f(to_sensor), 1.0f, std::numeric_limits<float>::infinity()};
      const std::uint8_t shadowed = scene.occluded(ray) ? 1 : 0;
      if (!(fl >= 0.0) || !(fb >= 0.0) || fl > static_cast<double>(lines - 1) ||
          fb > static_cast<double>(bins - 1)) {
        continue;
      }
      const auto l0 = std::min(static_cast<std::size_t>(fl), lines - 2);
      const auto b0 = std::min(static_cast<std::size_t>(fb), bins - 2);
      const double u = fb - static_cast<double>(b0);
      const double v = fl - static_cast<double>(l0);
      const auto& I = image.intensity;
      out.beta0.at(i, j) =
          static_cast<float>((1 - u) * (1 - v) * I.at(b0, l0) + u * (1 - v) * I.at(b0 + 1, l0) +
                             (1 - u) * v * I.at(b0, l0 + 1) + u * v * I.at(b0 + 1, l0 + 1));
      const auto ln = static_cast<std::size_t>(std::lround(fl));
      const auto bn = static_cast<std::size_t>(std::lround(fb));
      out.layover.at(i, j) = image.layover.at(bn, ln);
      out.shadow.at(i, j) = shadowed;
      out.valid.at(i, j) = 1;
    }
  });
  return out;
}

}  // namespace bsar
