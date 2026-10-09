#pragma once

#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

#include "backscatter/image/image.hpp"
#include "backscatter/scene/material.hpp"
#include "backscatter/scene/scene.hpp"
#include "backscatter/sensor/trajectory.hpp"

namespace bsar {

struct GeometricConfig {
  double range_spacing = 5.0;    // slant-range bin size [m]
  double azimuth_spacing = 5.0;  // along-track distance between lines [m]
  double rays_per_bin = 4.0;     // angular oversampling relative to one range bin
  int max_bounces = 3;
  Polarization polarization = Polarization::VV;
  unsigned threads = 0;  // 0 = all hardware threads
  /// Two consecutive first hits further apart than this many expected ray
  /// spacings are treated as an occlusion boundary (start of shadow).
  double occlusion_jump_factor = 8.0;
  /// Optional explicit windows; NaN selects them from the scene bounds.
  double t_start = std::numeric_limits<double>::quiet_NaN();
  double t_end = std::numeric_limits<double>::quiet_NaN();
  double near_range = std::numeric_limits<double>::quiet_NaN();
  double far_range = std::numeric_limits<double>::quiet_NaN();
};

/// Slant-range / azimuth image. Column x is range bin, row y is azimuth line.
struct GeometricImage {
  Image<float> intensity;            // radar brightness (beta0-like), all bounces
  std::vector<Image<float>> bounce;  // bounce[k]: contribution of (k+1)-bounce paths
  Image<std::uint8_t> layover;       // 1 where a bin receives returns from >1 surface part
  Image<std::uint8_t> shadow;        // 1 where a bin inside the swath receives no return
  double near_range = 0.0;           // slant range of the start of bin 0 [m]
  double range_spacing = 0.0;
  double t_start = 0.0;        // azimuth time of line 0 [s]
  double line_interval = 0.0;  // [s]
  std::size_t rays_traced = 0;
};

/// Geometric/radiometric SAR rendering. Rays are cast in each line's
/// zero-Doppler plane, followed through up to `max_bounces` specular
/// reflections, and every hit that sees the sensor is binned by its
/// equivalent slant range (half the round-trip path length).
GeometricImage render_geometric(const Scene& scene, const Trajectory& trajectory,
                                const GeometricConfig& config, const ScatteringModel& model);
GeometricImage render_geometric(const Scene& scene, const Trajectory& trajectory,
                                const GeometricConfig& config);

/// Coarse time search followed by Newton refinement: the zero-Doppler time of
/// `target` on `trajectory`, robust to poor initial guesses on long orbits.
double find_zero_doppler_time(const Trajectory& trajectory, const Vec3d& target);

}  // namespace bsar
