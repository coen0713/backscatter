#pragma once

#include <array>
#include <cstdint>

namespace bsar {

/// Philox4x32-10 counter-based generator (Salmon et al., "Parallel random
/// numbers: as easy as 1, 2, 3", SC 2011). Output is a pure function of
/// (counter, key), so random draws keyed by e.g. a facet ID are identical no
/// matter which thread produces them or in which order.
class Philox4x32 {
 public:
  using Counter = std::array<std::uint32_t, 4>;
  using Key = std::array<std::uint32_t, 2>;

  static constexpr Counter generate(Counter ctr, Key key) {
    for (int round = 0; round < 10; ++round) {
      if (round > 0) {
        key[0] += kW0;
        key[1] += kW1;
      }
      const std::uint64_t p0 = std::uint64_t{kM0} * ctr[0];
      const std::uint64_t p1 = std::uint64_t{kM1} * ctr[2];
      const auto hi0 = static_cast<std::uint32_t>(p0 >> 32);
      const auto lo0 = static_cast<std::uint32_t>(p0);
      const auto hi1 = static_cast<std::uint32_t>(p1 >> 32);
      const auto lo1 = static_cast<std::uint32_t>(p1);
      ctr = {hi1 ^ ctr[1] ^ key[0], lo1, hi0 ^ ctr[3] ^ key[1], lo0};
    }
    return ctr;
  }

 private:
  static constexpr std::uint32_t kM0 = 0xD2511F53u;
  static constexpr std::uint32_t kM1 = 0xCD9E8D57u;
  static constexpr std::uint32_t kW0 = 0x9E3779B9u;
  static constexpr std::uint32_t kW1 = 0xBB67AE85u;
};

/// Map 32 random bits to a double in the open interval (0, 1).
constexpr double to_unit_open(std::uint32_t bits) {
  return (static_cast<double>(bits) + 0.5) * (1.0 / 4294967296.0);
}

/// Convenience wrapper: a 64-bit seed is the key, and a pair of 64-bit
/// indices (e.g. facet ID, sample index) is the counter.
class CounterRng {
 public:
  constexpr explicit CounterRng(std::uint64_t seed)
      : key_{static_cast<std::uint32_t>(seed), static_cast<std::uint32_t>(seed >> 32)} {}

  /// Four independent uniforms in (0, 1) for the given (a, b) index pair.
  [[nodiscard]] constexpr std::array<double, 4> uniform4(std::uint64_t a, std::uint64_t b) const {
    const auto out =
        Philox4x32::generate({static_cast<std::uint32_t>(a), static_cast<std::uint32_t>(a >> 32),
                              static_cast<std::uint32_t>(b), static_cast<std::uint32_t>(b >> 32)},
                             key_);
    return {to_unit_open(out[0]), to_unit_open(out[1]), to_unit_open(out[2]), to_unit_open(out[3])};
  }

 private:
  Philox4x32::Key key_;
};

}  // namespace bsar
