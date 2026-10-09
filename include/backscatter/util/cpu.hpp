#pragma once

namespace bsar {

/// True when the CPU and OS support AVX2 and FMA (checked once, via CPUID and
/// XGETBV). SIMD kernels are compiled for these extensions per function and
/// only called when this returns true, so the default build still runs on
/// older x86-64 CPUs.
bool cpu_has_avx2_fma();

}  // namespace bsar

#if defined(__x86_64__) || defined(_M_X64)
#define BSAR_X86_64 1
#else
#define BSAR_X86_64 0
#endif

// Enable AVX2/FMA code generation for a single function. MSVC allows the
// intrinsics anywhere, so the attribute is only needed for GCC and Clang.
#if BSAR_X86_64 && (defined(__GNUC__) || defined(__clang__))
#define BSAR_TARGET_AVX2 __attribute__((target("avx2,fma")))
#else
#define BSAR_TARGET_AVX2
#endif
