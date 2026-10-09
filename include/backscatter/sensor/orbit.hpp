#pragma once

#include <filesystem>
#include <span>
#include <string_view>
#include <vector>

#include "backscatter/math/vec3.hpp"

namespace bsar {

/// One orbit state vector in an Earth-fixed frame. Time is in seconds since
/// 2000-01-01T00:00:00 UTC (see parse_utc).
struct StateVector {
  double t = 0.0;
  Vec3d position;
  Vec3d velocity;
};

/// Sorted state vectors with piecewise cubic Hermite interpolation. Using both
/// position and velocity at each node gives sub-millimetre position error for
/// Sentinel-1's 10 s vector spacing.
class Orbit {
 public:
  Orbit() = default;
  explicit Orbit(std::vector<StateVector> vectors);

  [[nodiscard]] Vec3d position(double t) const;
  [[nodiscard]] Vec3d velocity(double t) const;
  [[nodiscard]] double t_begin() const;
  [[nodiscard]] double t_end() const;
  [[nodiscard]] std::span<const StateVector> state_vectors() const { return vectors_; }

 private:
  [[nodiscard]] std::size_t segment(double t) const;
  std::vector<StateVector> vectors_;
};

/// Parse "YYYY-MM-DDThh:mm:ss[.ffffff]" (optionally prefixed "UTC=") into
/// seconds since 2000-01-01T00:00:00 UTC. Leap seconds are ignored, which is
/// consistent within one acquisition.
double parse_utc(std::string_view text);

/// Parse a Sentinel-1 orbit file (POEORB/RESORB, Earth Explorer XML) by
/// reading every <OSV> block.
Orbit parse_eof(std::string_view xml);
Orbit load_eof(const std::filesystem::path& path);

}  // namespace bsar
