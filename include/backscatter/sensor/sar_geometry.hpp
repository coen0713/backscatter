#pragma once

#include "backscatter/math/constants.hpp"
#include "backscatter/math/vec3.hpp"
#include "backscatter/sensor/trajectory.hpp"

namespace bsar {

struct ZeroDopplerSolution {
  double t = 0.0;      // azimuth time of closest approach
  double range = 0.0;  // one-way slant range at that time [m]
  bool converged = false;
  int iterations = 0;
};

/// Find the azimuth time at which the sensor velocity is perpendicular to the
/// line of sight to `target` (zero Doppler), by Newton iteration on
/// f(t) = v(t) . (p(t) - x). The derivative uses a central difference of the
/// velocity for the acceleration term.
ZeroDopplerSolution solve_zero_doppler(const Trajectory& trajectory, const Vec3d& target,
                                       double t_guess, double tolerance = 1e-9,
                                       int max_iterations = 50);

/// Theoretical slant-range resolution of an unweighted chirp, c / (2B).
/// (The -3 dB width of the sinc response is 0.886 times this.)
constexpr double slant_range_resolution(double bandwidth) {
  return kSpeedOfLight / (2.0 * bandwidth);
}

/// Stripmap azimuth resolution for an antenna of length L, L / 2.
constexpr double stripmap_azimuth_resolution(double antenna_length) { return antenna_length / 2.0; }

/// -3 dB width of a sinc-shaped impulse response relative to its nominal
/// resolution (1/bandwidth): 0.8859.
inline constexpr double kSinc3dbFactor = 0.885892941378904;

}  // namespace bsar
