#pragma once

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "backscatter/math/fft.hpp"
#include "backscatter/scene/material.hpp"
#include "backscatter/scene/scene.hpp"
#include "backscatter/sensor/radar.hpp"
#include "backscatter/sensor/trajectory.hpp"

namespace bsar {

/// A scattering path with `amplitude` = sqrt(radar cross-section) in metres.
/// The wave enters at `position` and leaves towards the sensor from `exit`
/// after travelling `internal_path` metres between bounces, so its one-way
/// equivalent range from sensor position p is
///   R(p) = (|p - position| + internal_path + |p - exit|) / 2.
/// A single-bounce point scatterer has exit == position and internal_path 0.
struct Scatterer {
  Vec3d position;
  double amplitude = 1.0;
  Vec3d exit;
  double internal_path = 0.0;
  int bounces = 1;

  Scatterer() = default;
  Scatterer(const Vec3d& p, double a) : position(p), amplitude(a), exit(p) {}
  Scatterer(const Vec3d& entry, const Vec3d& exit_point, double path, double a, int n)
      : position(entry), amplitude(a), exit(exit_point), internal_path(path), bounces(n) {}

  [[nodiscard]] double range_from(const Vec3d& p) const {
    return 0.5 * (length(p - position) + internal_path + length(p - exit));
  }
};

struct ScattererConfig {
  double density = 1.0;  // scatterers per square metre of surface
  std::uint64_t seed = 1;
  Polarization polarization = Polarization::VV;
  bool check_visibility = true;  // drop scatterers occluded from the illuminator
  int max_bounces = 3;           // 1 = single bounce only
};

/// Distribute sub-resolution scatterers over every facet. Positions are drawn
/// from a counter-based RNG keyed by (seed, facet id, sample index), so the
/// set is identical for any thread count or facet ordering. Each facet's
/// sigma0 (from the scattering model at its local incidence, seen from
/// `illuminator`) is split evenly across its scatterers; speckle then comes
/// from the random sub-wavelength positions alone.
///
/// With max_bounces > 1, the specular reflection of the illuminating ray at
/// each scatterer is traced on; every further hit that sees the illuminator
/// adds a multi-bounce path (entry, exit, internal length) weighted the same
/// way as in the geometric integrator. This is what produces the bright
/// double-bounce line at the foot of a wall. Bounce points are held fixed
/// over the aperture (they are found from `illuminator`), which is exact for
/// dihedral corners parallel to the track and an approximation otherwise.
std::vector<Scatterer> sample_scatterers(const Scene& scene, const Vec3d& illuminator,
                                         const ScattererConfig& config,
                                         const ScatteringModel& model);

/// Raw (uncompressed) echoes. Sample n of pulse k is at fast time
/// 2 * near_range / c + n / sample_rate.
struct RawData {
  std::size_t num_pulses = 0;
  std::size_t num_samples = 0;
  std::vector<cfloat> samples;  // num_pulses x num_samples
  std::vector<double> times;
  std::vector<Vec3d> positions;
  double near_range = 0.0;
  std::size_t window_samples = 0;  // samples spanning [near_range, far_range]
  RadarParams radar;
};

struct EchoWindow {
  double t_start = 0.0;
  double t_end = 0.0;
  double near_range = 0.0;
  double far_range = 0.0;
};

/// Sum delayed, phase-shifted chirps from every scatterer inside the azimuth
/// beam, per pulse. Each echo is placed with a band-limited (Kaiser-windowed
/// sinc) fractional delay and then convolved with the chirp by FFT. Pulses
/// are independent, so the result does not depend on the thread count.
RawData synthesize_raw(std::span<const Scatterer> scatterers, const Trajectory& trajectory,
                       const RadarParams& radar, const EchoWindow& window, unsigned threads = 0);

/// Range-compressed echoes: complex samples spaced `range_spacing` apart in
/// one-way slant range, starting at `near_range`.
struct CompressedData {
  std::size_t num_pulses = 0;
  std::size_t num_samples = 0;
  std::vector<cfloat> samples;
  std::vector<Vec3d> positions;
  double near_range = 0.0;
  double range_spacing = 0.0;
  double wavelength = 0.0;
};

/// Matched filtering with the transmitted chirp (FFT correlation), upsampled
/// by `upsample` through zero-padding in frequency so backprojection can use
/// linear interpolation. The output is normalised so a unit scatterer
/// compresses to a unit peak.
CompressedData range_compress(const RawData& raw, int upsample = 8, unsigned threads = 0);

/// Pick an echo window that covers the scene for its full synthetic aperture.
EchoWindow auto_echo_window(const Aabb& bounds, const Trajectory& trajectory,
                            const RadarParams& radar);

}  // namespace bsar
