#include "backscatter/sensor/sar_geometry.hpp"

#include <cmath>
#include <stdexcept>

namespace bsar {

ZeroDopplerSolution solve_zero_doppler(const Trajectory& trajectory, const Vec3d& target,
                                       double t_guess, double tolerance, int max_iterations) {
  ZeroDopplerSolution sol;
  double t = t_guess;
  constexpr double kDiffStep = 1e-3;  // s, for the acceleration estimate
  for (int i = 0; i < max_iterations; ++i) {
    const Vec3d p = trajectory.position(t);
    const Vec3d v = trajectory.velocity(t);
    const Vec3d a =
        (trajectory.velocity(t + kDiffStep) - trajectory.velocity(t - kDiffStep)) / (2 * kDiffStep);
    const Vec3d los = p - target;
    const double f = dot(v, los);
    const double df = dot(a, los) + dot(v, v);
    if (df == 0.0) {
      break;
    }
    const double step = f / df;
    t -= step;
    sol.iterations = i + 1;
    if (std::abs(step) < tolerance) {
      sol.converged = true;
      break;
    }
  }
  sol.t = t;
  sol.range = length(trajectory.position(t) - target);
  return sol;
}

LinearTrajectory make_side_looking_track(const Aabb& scene, double altitude, double speed,
                                         double incidence_deg, LookSide side) {
  if (scene.empty()) {
    throw std::invalid_argument("make_side_looking_track: empty scene bounds");
  }
  const Vec3d c(scene.center());
  const double ground_offset = altitude * std::tan(deg_to_rad(incidence_deg));
  // Heading north, a right-looking radar sees targets to the east, so it
  // flies west of the scene.
  const double x = side == LookSide::Right ? c.x - ground_offset : c.x + ground_offset;
  const Vec3d p0{x, c.y, c.z + altitude};
  const double half_span = 1e5;  // s; long enough for any scene
  return LinearTrajectory(p0, Vec3d{0.0, speed, 0.0}, -half_span, half_span);
}

}  // namespace bsar
