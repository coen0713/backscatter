#include "backscatter/util/cpu.hpp"

#if BSAR_X86_64
#if defined(_MSC_VER) && !defined(__clang__)
#include <immintrin.h>
#include <intrin.h>
#else
#include <cpuid.h>
#endif
#endif

namespace bsar {
namespace {

bool detect() {
#if BSAR_X86_64
  unsigned a = 0;
  unsigned b = 0;
  unsigned c = 0;
  unsigned d = 0;
#if defined(_MSC_VER) && !defined(__clang__)
  int info[4];
  __cpuid(info, 0);
  if (info[0] < 7) {
    return false;
  }
  __cpuid(info, 1);
  c = static_cast<unsigned>(info[2]);
#else
  if (__get_cpuid_max(0, nullptr) < 7 || __get_cpuid(1, &a, &b, &c, &d) == 0) {
    return false;
  }
#endif
  const bool fma = (c & (1u << 12)) != 0;
  const bool osxsave = (c & (1u << 27)) != 0;
  const bool avx = (c & (1u << 28)) != 0;
  if (!fma || !osxsave || !avx) {
    return false;
  }
  // The OS must save the YMM registers on context switches (XCR0 bits 1, 2).
#if defined(_MSC_VER) && !defined(__clang__)
  const unsigned long long xcr0 = _xgetbv(0);
#else
  unsigned lo = 0;
  unsigned hi = 0;
  __asm__ volatile("xgetbv" : "=a"(lo), "=d"(hi) : "c"(0));
  const unsigned long long xcr0 = (static_cast<unsigned long long>(hi) << 32) | lo;
#endif
  if ((xcr0 & 6u) != 6u) {
    return false;
  }
#if defined(_MSC_VER) && !defined(__clang__)
  __cpuidex(info, 7, 0);
  b = static_cast<unsigned>(info[1]);
#else
  __cpuid_count(7, 0, a, b, c, d);
#endif
  return (b & (1u << 5)) != 0;  // AVX2
#else
  return false;
#endif
}

}  // namespace

bool cpu_has_avx2_fma() {
  static const bool has = detect();
  return has;
}

}  // namespace bsar
