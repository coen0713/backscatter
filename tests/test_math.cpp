#include <atomic>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <vector>

#include "backscatter/math/constants.hpp"
#include "backscatter/math/fft.hpp"
#include "backscatter/math/rng.hpp"
#include "backscatter/math/vec3.hpp"
#include "backscatter/util/parallel.hpp"

using namespace bsar;
using Catch::Approx;

TEST_CASE("Vec3 basic algebra", "[math]") {
  const Vec3d a{1, 2, 3};
  const Vec3d b{4, -5, 6};
  CHECK(dot(a, b) == Approx(12.0));
  const Vec3d c = cross(a, b);
  CHECK(dot(c, a) == Approx(0.0).margin(1e-12));
  CHECK(dot(c, b) == Approx(0.0).margin(1e-12));
  CHECK(length(normalize(b)) == Approx(1.0));
  const Vec3d n{0, 0, 1};
  const Vec3d r = reflect(Vec3d{1, 0, -1}, n);
  CHECK(r == Vec3d{1, 0, 1});
  CHECK(max_dimension(Vec3d{-1, 5, 2}) == 1);
}

TEST_CASE("Philox4x32-10 matches the Random123 known-answer vectors", "[math][rng]") {
  const auto zero = Philox4x32::generate({0, 0, 0, 0}, {0, 0});
  CHECK(zero == Philox4x32::Counter{0x6627e8d5u, 0xe169c58du, 0xbc57ac4cu, 0x9b00dbd8u});
  const auto ones = Philox4x32::generate({0xffffffffu, 0xffffffffu, 0xffffffffu, 0xffffffffu},
                                         {0xffffffffu, 0xffffffffu});
  CHECK(ones == Philox4x32::Counter{0x408f276du, 0x41c83b0eu, 0xa20bc7c6u, 0x6d5451fdu});
  const auto pi = Philox4x32::generate({0x243f6a88u, 0x85a308d3u, 0x13198a2eu, 0x03707344u},
                                       {0xa4093822u, 0x299f31d0u});
  CHECK(pi == Philox4x32::Counter{0xd16cfe09u, 0x94fdccebu, 0x5001e420u, 0x24126ea1u});
}

TEST_CASE("CounterRng is uniform and keyed", "[math][rng]") {
  const CounterRng rng(42);
  double sum = 0.0;
  constexpr int kN = 20000;
  for (int i = 0; i < kN; ++i) {
    for (const double u : rng.uniform4(static_cast<std::uint64_t>(i), 7)) {
      REQUIRE(u > 0.0);
      REQUIRE(u < 1.0);
      sum += u;
    }
  }
  CHECK(sum / (4.0 * kN) == Approx(0.5).margin(0.005));
  CHECK(rng.uniform4(3, 9) == rng.uniform4(3, 9));
  CHECK(rng.uniform4(3, 9) != CounterRng(43).uniform4(3, 9));
}

TEST_CASE("FFT agrees with a direct DFT and inverts", "[math][fft]") {
  for (const std::size_t n : {1u, 2u, 8u, 64u, 1024u}) {
    std::vector<cdouble> x(n);
    for (std::size_t i = 0; i < n; ++i) {
      x[i] = {std::sin(0.3 * static_cast<double>(i)) + 0.1 * static_cast<double>(i % 7),
              std::cos(1.7 * static_cast<double>(i))};
    }
    std::vector<cdouble> X = x;
    fft_inplace(X);
    if (n <= 64) {
      for (std::size_t k = 0; k < n; ++k) {
        cdouble ref{};
        for (std::size_t j = 0; j < n; ++j) {
          ref += x[j] *
                 std::polar(1.0, -2.0 * kPi * static_cast<double>(j * k) / static_cast<double>(n));
        }
        REQUIRE(std::abs(X[k] - ref) < 1e-9);
      }
    }
    fft_inplace(X, true);
    for (std::size_t i = 0; i < n; ++i) {
      REQUIRE(std::abs(X[i] - x[i]) < 1e-12);
    }
  }
  std::vector<cdouble> bad(3);
  CHECK_THROWS(fft_inplace(bad));
  CHECK(next_pow2(1000) == 1024);
  CHECK(next_pow2(1024) == 1024);
}

TEST_CASE("parallel_for visits every index exactly once", "[util]") {
  std::vector<std::atomic<int>> hits(10007);
  parallel_for(0, hits.size(), 8, [&](std::size_t i) { hits[i].fetch_add(1); }, 13);
  for (const auto& h : hits) {
    REQUIRE(h.load() == 1);
  }
  CHECK_THROWS(parallel_for(0, 100, 4, [](std::size_t i) {
    if (i == 50) {
      throw std::runtime_error("boom");
    }
  }));
}
