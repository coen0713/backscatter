#pragma once

#include <complex>
#include <cstddef>
#include <span>

namespace bsar {

using cdouble = std::complex<double>;
using cfloat = std::complex<float>;

[[nodiscard]] constexpr bool is_pow2(std::size_t n) { return n != 0 && (n & (n - 1)) == 0; }

[[nodiscard]] constexpr std::size_t next_pow2(std::size_t n) {
  std::size_t p = 1;
  while (p < n) {
    p <<= 1;
  }
  return p;
}

/// In-place iterative radix-2 FFT. `data.size()` must be a power of two.
/// The forward transform uses exp(-2*pi*i*k*n/N); the inverse is scaled by
/// 1/N so that ifft(fft(x)) == x.
void fft_inplace(std::span<cdouble> data, bool inverse = false);

}  // namespace bsar
