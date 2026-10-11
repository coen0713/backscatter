#pragma once

#include <string>
#include <vector>

#include "backscatter/sensor/geodesy.hpp"
#include "backscatter/sensor/orbit.hpp"
#include "backscatter/sensor/trajectory.hpp"

namespace bsar {

/// One closest approach of an orbit to a ground point.
struct OrbitPass {
  double t = 0.0;              // zero-Doppler time, seconds since 2000-01-01 UTC
  double range = 0.0;          // slant range at zero Doppler [m]
  double incidence_deg = 0.0;  // at the target, from the local vertical
  bool ascending = false;      // satellite moving north (ECEF velocity z > 0)
  LookSide side = LookSide::Right;
};

/// Every closest approach of `orbit` to `target` within `max_range`, in time
/// order. Sentinel-1 is right-looking with IW incidence of roughly 30-46 deg,
/// so those are the passes that correspond to real acquisitions.
std::vector<OrbitPass> find_passes(const Orbit& orbit, const Geodetic& target,
                                   double max_range = 1.5e6);

/// Inverse of parse_utc: "YYYY-MM-DDThh:mm:ss.ffffff".
std::string format_utc(double seconds_since_2000);

}  // namespace bsar
