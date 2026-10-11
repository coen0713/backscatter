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
///
/// Time precision: absolute times (~7.7e8 s since 2000) only resolve ~1.2e-7 s
/// in a double, i.e. ~1 mm of satellite motion, which is ~0.2 rad of two-way
/// C-band phase. The orbit therefore works in *orbit time*: seconds since its
/// first state vector (`epoch()`, absolute). position(), velocity(), t_begin()
/// and t_end() all use orbit time; convert with `epoch()` at the edges only.
class Orbit {
 public:
  Orbit() = default;
  /// `vectors` carry absolute times (seconds since 2000-01-01 UTC).
  explicit Orbit(std::vector<StateVector> vectors);

  [[nodiscard]] Vec3d position(double t) const;
  [[nodiscard]] Vec3d velocity(double t) const;
  [[nodiscard]] double t_begin() const;  // orbit time, always 0
  [[nodiscard]] double t_end() const;    // orbit time
  /// Absolute time of orbit time 0 (seconds since 2000-01-01 UTC).
  [[nodiscard]] double epoch() const { return epoch_; }
  /// State vectors with times in orbit time.
  [[nodiscard]] std::span<const StateVector> state_vectors() const { return vectors_; }
  /// A copy restricted to state vectors within `half_width` seconds of the
  /// absolute time `centre` (keeps the same kind of interpolation).
  [[nodiscard]] Orbit window(double centre, double half_width) const;

 private:
  [[nodiscard]] std::size_t segment(double t) const;
  std::vector<StateVector> vectors_;
  double epoch_ = 0.0;
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
