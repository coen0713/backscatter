#include "backscatter/math/fft.hpp"

#include <stdexcept>
#include <utility>
#include <vector>

#include "backscatter/math/constants.hpp"

namespace bsar {

void fft_inplace(std::span<cdouble> data, bool inverse) {
  const std::size_t n = data.size();
  if (n <= 1) {
    return;
  }
  if (!is_pow2(n)) {
    throw std::invalid_argument("fft_inplace: size must be a power of two");
  }

  // Bit-reversal permutation.
  for (std::size_t i = 1, j = 0; i < n; ++i) {
    std::size_t bit = n >> 1;
    for (; (j & bit) != 0; bit >>= 1) {
      j ^= bit;
    }
    j ^= bit;
    if (i < j) {
      std::swap(data[i], data[j]);
    }
  }

  // Twiddles computed directly (not by repeated multiplication) to keep the
  // error at O(eps * log N) for long transforms.
  const double sign = inverse ? 1.0 : -1.0;
  std::vector<cdouble> twiddle(n / 2);
  for (std::size_t k = 0; k < n / 2; ++k) {
    twiddle[k] =
        std::polar(1.0, sign * 2.0 * kPi * static_cast<double>(k) / static_cast<double>(n));
  }

  for (std::size_t len = 2; len <= n; len <<= 1) {
    const std::size_t half = len / 2;
    const std::size_t stride = n / len;
    for (std::size_t start = 0; start < n; start += len) {
      for (std::size_t k = 0; k < half; ++k) {
        const cdouble u = data[start + k];
        const cdouble v = data[start + k + half] * twiddle[k * stride];
        data[start + k] = u + v;
        data[start + k + half] = u - v;
      }
    }
  }

  if (inverse) {
    const double scale = 1.0 / static_cast<double>(n);
    for (auto& x : data) {
      x *= scale;
    }
  }
}

}  // namespace bsar
