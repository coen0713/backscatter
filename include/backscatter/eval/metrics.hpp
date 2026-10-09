#pragma once

#include <cstddef>
#include <cstdint>
#include <span>

namespace bsar {

/// Impulse response metrics from a 1-D amplitude cut through a point target.
struct IrfMetrics {
  double peak_position = 0.0;  // in samples (sub-sample, parabolic fit)
  double peak_amplitude = 0.0;
  double resolution_3db = 0.0;  // -3 dB (half-power) width, in units of `spacing`
  double pslr_db = 0.0;         // peak sidelobe ratio
  double islr_db = 0.0;         // integrated sidelobe ratio (sidelobes / main lobe)
  bool valid = false;
};

/// The main lobe extends from the peak out to the first minimum on each side.
IrfMetrics analyze_irf(std::span<const float> amplitude, double spacing);

/// Single-look speckle statistics for intensity samples.
struct SpeckleStats {
  std::size_t n = 0;
  double mean = 0.0;
  double variance = 0.0;
  double enl = 0.0;  // equivalent number of looks, mean^2 / variance
  /// Kolmogorov-Smirnov distance to an exponential with the sample mean, and
  /// Stephens' (1974) modified statistic for the estimated-mean case.
  double ks_d = 0.0;
  double ks_modified = 0.0;
  bool exponential_at_5pct = false;  // ks_modified < 1.094
  bool exponential_at_1pct = false;  // ks_modified < 1.308
};

SpeckleStats speckle_statistics(std::span<const double> intensity);

/// Intersection over union of two binary masks (non-zero = set). Returns 1
/// when both masks are empty.
double mask_iou(std::span<const std::uint8_t> a, std::span<const std::uint8_t> b);

}  // namespace bsar
