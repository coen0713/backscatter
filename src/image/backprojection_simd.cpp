// AVX2/FMA backprojection: four pixels per step, double-precision range and
// phase, polynomial sincos. Kept free of standard-library templates so no
// AVX2-encoded copy of a shared inline function can leak into other code.

#include "backprojection_simd.hpp"

#include "backscatter/util/cpu.hpp"

#if BSAR_X86_64
#include <immintrin.h>
#endif

namespace bsar::detail {

#if BSAR_X86_64

namespace {

/// sin and cos of `phase` (any magnitude up to ~1e9 rad). The argument is
/// reduced to [-pi/4, pi/4] with a two-part pi/2 and an exact FMA, then
/// Taylor polynomials (error < 5e-12) are combined by quadrant.
BSAR_TARGET_AVX2 inline void sincos4(__m256d phase, __m256d& s_out, __m256d& c_out) {
  const __m256d two_over_pi = _mm256_set1_pd(0.63661977236758134308);
  const __m256d pio2_hi = _mm256_set1_pd(1.5707963267948966192);
  const __m256d pio2_lo = _mm256_set1_pd(6.123233995736766036e-17);
  const __m256d q = _mm256_round_pd(_mm256_mul_pd(phase, two_over_pi),
                                    _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
  __m256d x = _mm256_fnmadd_pd(q, pio2_hi, phase);
  x = _mm256_fnmadd_pd(q, pio2_lo, x);
  const __m256d x2 = _mm256_mul_pd(x, x);

  __m256d ps = _mm256_set1_pd(-1.0 / 39916800.0);
  ps = _mm256_fmadd_pd(ps, x2, _mm256_set1_pd(1.0 / 362880.0));
  ps = _mm256_fmadd_pd(ps, x2, _mm256_set1_pd(-1.0 / 5040.0));
  ps = _mm256_fmadd_pd(ps, x2, _mm256_set1_pd(1.0 / 120.0));
  ps = _mm256_fmadd_pd(ps, x2, _mm256_set1_pd(-1.0 / 6.0));
  const __m256d sin_x = _mm256_fmadd_pd(_mm256_mul_pd(ps, x2), x, x);

  __m256d pc = _mm256_set1_pd(1.0 / 479001600.0);
  pc = _mm256_fmadd_pd(pc, x2, _mm256_set1_pd(-1.0 / 3628800.0));
  pc = _mm256_fmadd_pd(pc, x2, _mm256_set1_pd(1.0 / 40320.0));
  pc = _mm256_fmadd_pd(pc, x2, _mm256_set1_pd(-1.0 / 720.0));
  pc = _mm256_fmadd_pd(pc, x2, _mm256_set1_pd(1.0 / 24.0));
  pc = _mm256_fmadd_pd(pc, x2, _mm256_set1_pd(-0.5));
  const __m256d cos_x = _mm256_fmadd_pd(pc, x2, _mm256_set1_pd(1.0));

  // Quadrant q mod 4: odd quadrants swap sin/cos; signs follow
  // sin: q & 2, cos: (q + 1) & 2.
  const __m128i qi = _mm256_cvtpd_epi32(q);
  const __m128i one = _mm_set1_epi32(1);
  const __m128i two = _mm_set1_epi32(2);
  const __m256d swap =
      _mm256_castsi256_pd(_mm256_cvtepi32_epi64(_mm_cmpeq_epi32(_mm_and_si128(qi, one), one)));
  const __m256d sin_neg =
      _mm256_castsi256_pd(_mm256_cvtepi32_epi64(_mm_cmpeq_epi32(_mm_and_si128(qi, two), two)));
  const __m256d cos_neg = _mm256_castsi256_pd(
      _mm256_cvtepi32_epi64(_mm_cmpeq_epi32(_mm_and_si128(_mm_add_epi32(qi, one), two), two)));
  const __m256d sign_bit = _mm256_set1_pd(-0.0);
  const __m256d s = _mm256_blendv_pd(sin_x, cos_x, swap);
  const __m256d c = _mm256_blendv_pd(cos_x, sin_x, swap);
  s_out = _mm256_xor_pd(s, _mm256_and_pd(sin_neg, sign_bit));
  c_out = _mm256_xor_pd(c, _mm256_and_pd(cos_neg, sign_bit));
}

}  // namespace

BSAR_TARGET_AVX2 std::size_t accumulate_avx2(const BackprojectionInputs& in, const double* px,
                                             const double* py, const double* pz, std::size_t n,
                                             std::size_t k0, std::size_t k1, double* acc) {
  const std::size_t n4 = n & ~static_cast<std::size_t>(3);
  const __m256d near = _mm256_set1_pd(in.near_range);
  const __m256d inv = _mm256_set1_pd(in.inv_spacing);
  const __m256d phase_scale = _mm256_set1_pd(in.phase_scale);
  const __m256d last = _mm256_set1_pd(static_cast<double>(in.num_samples - 1));
  const __m256d zero = _mm256_setzero_pd();

  for (std::size_t i = 0; i < n4; i += 4) {
    const __m256d x = _mm256_loadu_pd(px + i);
    const __m256d y = _mm256_loadu_pd(py + i);
    const __m256d z = _mm256_loadu_pd(pz + i);
    __m256d acc01 = _mm256_loadu_pd(acc + 2 * i);      // re0 im0 re1 im1
    __m256d acc23 = _mm256_loadu_pd(acc + 2 * i + 4);  // re2 im2 re3 im3

    for (std::size_t k = k0; k < k1; ++k) {
      const double* p = in.pulse_xyz + 3 * k;
      const __m256d dx = _mm256_sub_pd(x, _mm256_set1_pd(p[0]));
      const __m256d dy = _mm256_sub_pd(y, _mm256_set1_pd(p[1]));
      const __m256d dz = _mm256_sub_pd(z, _mm256_set1_pd(p[2]));
      const __m256d r = _mm256_sqrt_pd(_mm256_add_pd(
          _mm256_add_pd(_mm256_mul_pd(dx, dx), _mm256_mul_pd(dy, dy)), _mm256_mul_pd(dz, dz)));
      const __m256d s = _mm256_mul_pd(_mm256_sub_pd(r, near), inv);
      const __m256d valid =
          _mm256_and_pd(_mm256_cmp_pd(s, zero, _CMP_GE_OQ), _mm256_cmp_pd(s, last, _CMP_LT_OQ));
      if (_mm256_movemask_pd(valid) == 0) {
        continue;
      }
      const __m256d sf = _mm256_floor_pd(s);
      const __m256d frac = _mm256_sub_pd(s, sf);
      const __m128i idx = _mm256_cvttpd_epi32(_mm256_blendv_pd(zero, sf, valid));

      // One complex<float> is 8 bytes, so gather it as a double.
      const auto* row = reinterpret_cast<const double*>(in.samples + k * in.num_samples);
      const __m256 v0 = _mm256_castpd_ps(_mm256_i32gather_pd(row, idx, 8));
      const __m256 v1 = _mm256_castpd_ps(_mm256_i32gather_pd(row + 1, idx, 8));
      const __m128 f4 = _mm256_cvtpd_ps(frac);
      const __m256 f8 = _mm256_set_m128(_mm_unpackhi_ps(f4, f4), _mm_unpacklo_ps(f4, f4));
      __m256 v = _mm256_fmadd_ps(_mm256_sub_ps(v1, v0), f8, v0);
      v = _mm256_and_ps(v, _mm256_castpd_ps(valid));
      const __m256d v01 = _mm256_cvtps_pd(_mm256_castps256_ps128(v));
      const __m256d v23 = _mm256_cvtps_pd(_mm256_extractf128_ps(v, 1));

      __m256d sn;
      __m256d cs;
      sincos4(_mm256_mul_pd(r, phase_scale), sn, cs);
      // (a + ib)(c + is) = (ac - bs) + i(as + bc) with interleaved (a, b).
      const __m256d c01 = _mm256_permute4x64_pd(cs, 0x50);  // c0 c0 c1 c1
      const __m256d c23 = _mm256_permute4x64_pd(cs, 0xFA);  // c2 c2 c3 c3
      const __m256d s01 = _mm256_permute4x64_pd(sn, 0x50);
      const __m256d s23 = _mm256_permute4x64_pd(sn, 0xFA);
      acc01 =
          _mm256_add_pd(acc01, _mm256_addsub_pd(_mm256_mul_pd(v01, c01),
                                                _mm256_mul_pd(_mm256_permute_pd(v01, 0x5), s01)));
      acc23 =
          _mm256_add_pd(acc23, _mm256_addsub_pd(_mm256_mul_pd(v23, c23),
                                                _mm256_mul_pd(_mm256_permute_pd(v23, 0x5), s23)));
    }
    _mm256_storeu_pd(acc + 2 * i, acc01);
    _mm256_storeu_pd(acc + 2 * i + 4, acc23);
  }
  return n4;
}

#else

std::size_t accumulate_avx2(const BackprojectionInputs&, const double*, const double*,
                            const double*, std::size_t, std::size_t, std::size_t, double*) {
  return 0;
}

#endif

}  // namespace bsar::detail
