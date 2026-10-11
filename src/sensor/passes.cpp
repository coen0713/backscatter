#include "backscatter/sensor/passes.hpp"

#include <cmath>
#include <cstdio>

#include "backscatter/math/constants.hpp"
#include "backscatter/sensor/sar_geometry.hpp"

namespace bsar {

std::vector<OrbitPass> find_passes(const Orbit& orbit, const Geodetic& target, double max_range) {
  // In the target's ENU frame the target is the origin, so the squared range
  // |p|^2 has derivative 2 v.p, which crosses zero from below at each closest
  // approach.
  const EnuFrame frame(target);
  const OrbitTrajectory traj(orbit, frame);
  std::vector<OrbitPass> passes;
  constexpr double kStep = 10.0;  // s; S1 passes last minutes, so no minimum is missed
  double t_prev = traj.t_begin();
  double f_prev = dot(traj.velocity(t_prev), traj.position(t_prev));
  for (double t = t_prev + kStep; t <= traj.t_end(); t += kStep) {
    const double f = dot(traj.velocity(t), traj.position(t));
    if (f_prev < 0.0 && f >= 0.0) {
      const ZeroDopplerSolution zd = solve_zero_doppler(traj, Vec3d{}, t - 0.5 * kStep);
      const Vec3d p = traj.position(zd.t);
      const Vec3d v = traj.velocity(zd.t);
      if (zd.converged && zd.range <= max_range) {
        OrbitPass pass;
        pass.t = orbit.epoch() + zd.t;  // absolute, for display and selection
        pass.range = zd.range;
        pass.incidence_deg = rad_to_deg(std::acos(p.z / length(p)));
        pass.ascending = orbit.velocity(zd.t).z > 0.0;
        // Target (origin) relative to the sensor is -p; it is on the right of
        // the ground track if it lies along v x up.
        pass.side =
            dot(cross(v, Vec3d{0.0, 0.0, 1.0}), -p) > 0.0 ? LookSide::Right : LookSide::Left;
        passes.push_back(pass);
      }
    }
    t_prev = t;
    f_prev = f;
  }
  return passes;
}

std::string format_utc(double s) {
  const double day_seconds = 86400.0;
  auto days = static_cast<long long>(std::floor(s / day_seconds));
  double rem = s - static_cast<double>(days) * day_seconds;
  // Days since 2000-01-01 -> civil date (Howard Hinnant's civil_from_days,
  // shifted from the 1970 epoch).
  long long z = days + 10957 + 719468;
  const long long era = (z >= 0 ? z : z - 146096) / 146097;
  const auto doe = static_cast<unsigned>(z - era * 146097);
  const unsigned yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
  const long long y = static_cast<long long>(yoe) + era * 400;
  const unsigned doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
  const unsigned mp = (5 * doy + 2) / 153;
  const unsigned d = doy - (153 * mp + 2) / 5 + 1;
  const unsigned m = mp < 10 ? mp + 3 : mp - 9;
  const long long year = y + (m <= 2 ? 1 : 0);
  const int hour = static_cast<int>(rem / 3600.0);
  rem -= hour * 3600.0;
  const int minute = static_cast<int>(rem / 60.0);
  rem -= minute * 60.0;
  char buf[48];
  std::snprintf(buf, sizeof buf, "%04lld-%02u-%02uT%02d:%02d:%09.6f", year, m, d, hour, minute,
                rem);
  return buf;
}

}  // namespace bsar
