#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <vector>

#include "backscatter/eval/metrics.hpp"
#include "backscatter/math/constants.hpp"
#include "backscatter/math/rng.hpp"
#include "backscatter/sensor/sar_geometry.hpp"

using namespace bsar;
using Catch::Approx;

TEST_CASE("IRF metrics of an ideal sinc", "[eval]") {
  constexpr double kSpacing = 0.05;  // 20 samples per resolution cell
  std::vector<float> amp;
  for (int i = -400; i <= 400; ++i) {
    const double x = i * kSpacing;
    amp.push_back(static_cast<float>(std::abs(x == 0.0 ? 1.0 : std::sin(kPi * x) / (kPi * x))));
  }
  const IrfMetrics m = analyze_irf(amp, kSpacing);
  REQUIRE(m.valid);
  CHECK(m.resolution_3db == Approx(kSinc3dbFactor).epsilon(0.002));
  CHECK(m.pslr_db == Approx(-13.26).margin(0.05));
  CHECK(m.islr_db == Approx(-9.9).margin(0.3));  // -9.7 dB untruncated; cut at +/- 20 cells
  CHECK(m.peak_position == Approx(400.0).margin(1e-6));
}

TEST_CASE("Speckle statistics accept exponential and reject Rayleigh amplitude", "[eval]") {
  const CounterRng rng(2024);
  std::vector<double> intensity;
  std::vector<double> amplitude;
  for (std::uint64_t i = 0; i < 2000; ++i) {
    const double v = -3.0 * std::log(rng.uniform4(i, 0)[0]);  // exponential, mean 3
    intensity.push_back(v);
    amplitude.push_back(std::sqrt(v));
  }
  const SpeckleStats s = speckle_statistics(intensity);
  CHECK(s.mean == Approx(3.0).epsilon(0.06));
  CHECK(s.enl == Approx(1.0).margin(0.1));
  CHECK(s.exponential_at_5pct);
  const SpeckleStats a = speckle_statistics(amplitude);
  CHECK(a.enl > 2.0);  // Rayleigh amplitude has ENL ~ 3.7
  CHECK_FALSE(a.exponential_at_1pct);
}

TEST_CASE("Mask IoU", "[eval]") {
  const std::vector<std::uint8_t> a{1, 1, 0, 0, 1};
  const std::vector<std::uint8_t> b{1, 0, 0, 1, 1};
  CHECK(mask_iou(a, b) == Approx(0.5));
  const std::vector<std::uint8_t> z(5, 0);
  CHECK(mask_iou(z, z) == 1.0);
  CHECK_THROWS(mask_iou(a, std::vector<std::uint8_t>(3)));
}
