#pragma once

// Internal interface of the AVX2 backprojection kernel.

#include <complex>
#include <cstddef>

namespace bsar::detail {

struct BackprojectionInputs {
  const std::complex<float>* samples = nullptr;  // num_pulses x num_samples
  std::size_t num_samples = 0;
  const double* pulse_xyz = nullptr;  // 3 doubles per pulse
  double near_range = 0.0;
  double inv_spacing = 0.0;  // 1 / range spacing
  double phase_scale = 0.0;  // 4 pi / lambda
};

/// Add pulses [k0, k1) to the first (n & ~3) pixels. Pixel coordinates are
/// separate x/y/z arrays; `acc` holds interleaved (re, im) doubles. Returns
/// the number of pixels handled; the caller does the remaining (< 4) pixels.
/// Requires cpu_has_avx2_fma().
std::size_t accumulate_avx2(const BackprojectionInputs& in, const double* px, const double* py,
                            const double* pz, std::size_t n, std::size_t k0, std::size_t k1,
                            double* acc);

}  // namespace bsar::detail
