#include <algorithm>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>

#include "backscatter/eval/theory_checks.hpp"
#include "backscatter/image/backprojection.hpp"
#include "backscatter/math/constants.hpp"
#include "backscatter/scene/scene.hpp"
#include "backscatter/sensor/sar_geometry.hpp"
#include "backscatter/trace/coherent.hpp"

using namespace bsar;
using Catch::Approx;

TEST_CASE("Range compression puts a scatterer at its range with unit peak", "[coherent]") {
  RadarParams radar = theory_radar();
  const LinearTrajectory track({0.0, 0.0, 0.0}, {0.0, 100.0, 0.0}, -1.0, 1.0);
  const Vec3d target{1234.567, 0.0, 0.0};
  const std::vector<Scatterer> s{{target, 2.0}};
  EchoWindow w{0.0, 0.0, 1200.0, 1300.0};
  const RawData raw = synthesize_raw(s, track, radar, w, 1);
  REQUIRE(raw.num_pulses == 1);
  const CompressedData comp = range_compress(raw, 8, 1);
  std::size_t peak = 0;
  for (std::size_t i = 1; i < comp.num_samples; ++i) {
    if (std::abs(comp.samples[i]) > std::abs(comp.samples[peak])) {
      peak = i;
    }
  }
  const double peak_range = comp.near_range + static_cast<double>(peak) * comp.range_spacing;
  CHECK(std::abs(peak_range - target.x) <= comp.range_spacing);
  CHECK(std::abs(comp.samples[peak]) == Approx(2.0).epsilon(0.01));
  // Phase at the peak is -4 pi R / lambda (mod 2 pi).
  const double expected_phase = std::remainder(-4.0 * kPi * target.x / radar.wavelength(), 2 * kPi);
  const double phase_err = std::remainder(std::arg(comp.samples[peak]) - expected_phase, 2.0 * kPi);
  // The compressed baseband pulse is a real sinc, so its phase is flat across
  // the main lobe; only interpolation ripple remains.
  CHECK(std::abs(phase_err) < 0.05);
}

TEST_CASE("Point target: resolution and sidelobes match theory", "[coherent][theory]") {
  const PointTargetReport rep = run_point_target_check(theory_radar(), CoherentGeometry{});
  INFO("range 3dB " << rep.range.resolution_3db << " expected " << rep.expected_range_3db);
  INFO("azimuth 3dB " << rep.azimuth.resolution_3db << " expected " << rep.expected_azimuth_3db);
  REQUIRE(rep.range.valid);
  REQUIRE(rep.azimuth.valid);
  CHECK(std::abs(rep.range_error_pct) < 5.0);
  CHECK(std::abs(rep.azimuth_error_pct) < 5.0);
  CHECK(rep.range.pslr_db == Approx(-13.26).margin(1.0));
  CHECK(rep.azimuth.pslr_db == Approx(-13.26).margin(1.0));
  CHECK(rep.range.islr_db < -8.0);
  // Focused at the true position.
  CHECK(std::abs(rep.range.peak_position - 480.0) < 1.0);
  CHECK(std::abs(rep.azimuth.peak_position - 480.0) < 1.0);
}

TEST_CASE("Speckle: single-look intensity is exponential (ENL ~ 1)", "[coherent][theory]") {
  const SpeckleReport rep = run_speckle_check(theory_radar(), CoherentGeometry{});
  INFO("n " << rep.stats.n << " ENL " << rep.stats.enl << " KS* " << rep.stats.ks_modified);
  CHECK(rep.stats.n > 300);
  CHECK(rep.stats.enl == Approx(1.0).margin(0.15));
  CHECK(rep.stats.exponential_at_1pct);
}

TEST_CASE("Scatterer sampling is deterministic and density-correct", "[coherent]") {
  const Scene scene = make_plane_scene(50.0);
  ScattererConfig cfg;
  cfg.density = 3.0;
  cfg.seed = 99;
  const FacetScatteringModel model;
  const Vec3d sensor{-4000.0, 0.0, 3000.0};
  const auto a = sample_scatterers(scene, sensor, cfg, model);
  const auto b = sample_scatterers(scene, sensor, cfg, model);
  REQUIRE(a.size() == b.size());
  CHECK(static_cast<double>(a.size()) == Approx(3.0 * 2500.0).epsilon(0.01));
  for (std::size_t i = 0; i < a.size(); ++i) {
    REQUIRE(a[i].position == b[i].position);
  }
  // Seen from below the plane, nothing is visible.
  CHECK(sample_scatterers(scene, {0.0, 0.0, -1000.0}, cfg, model).empty());
  // Total RCS matches sigma0 x area.
  double rcs = 0.0;
  for (const auto& s : a) {
    rcs += s.amplitude * s.amplitude;
  }
  const double cos_inc = 3000.0 / std::hypot(4000.0, 3000.0);
  CHECK(rcs == Approx(model.backscatter(materials::dry_soil(), cos_inc, Polarization::VV) * 2500.0)
                   .epsilon(0.02));
}

TEST_CASE("Coherent dihedral focuses a double-bounce line at the wall foot", "[coherent]") {
  const RadarParams radar = theory_radar();
  const Scene scene = make_dihedral_scene(60.0, 10.0, 10.0, 8.0);  // west wall at x = -5
  const LinearTrajectory track({-4000.0, 0.0, 3000.0}, {0.0, 100.0, 0.0}, -100.0, 100.0);
  ScattererConfig sc;
  sc.density = 4.0;
  const auto scatterers = sample_scatterers(scene, track.position(0.0), sc, FacetScatteringModel{});
  const auto doubles = std::count_if(scatterers.begin(), scatterers.end(),
                                     [](const Scatterer& s) { return s.bounces == 2; });
  CHECK(doubles > 100);
  const EchoWindow w = auto_echo_window(scene.bounds(), track, radar);
  const CompressedData comp = range_compress(synthesize_raw(scatterers, track, radar, w), 8);

  // Ground-range cut through the building's middle, west of the wall.
  PixelGrid grid;
  grid.width = 241;
  grid.height = 1;
  grid.origin = {-29.0, 0.0, 0.0};
  grid.step_x = {0.1, 0.0, 0.0};
  const auto img = backproject(comp, grid);
  std::vector<double> power(img.size());
  for (std::size_t i = 0; i < img.size(); ++i) {
    power[i] = std::norm(std::complex<double>(img.data[i]));
  }
  const auto peak = std::max_element(power.begin(), power.end()) - power.begin();
  const double peak_x = -29.0 + 0.1 * static_cast<double>(peak);
  INFO("peak at x = " << peak_x);
  CHECK(std::abs(peak_x - (-5.0)) < 1.0);
  // Far above the open-ground speckle (x < -20 is beyond the wall's layover).
  double ground = 0.0;
  for (std::size_t i = 0; i < 60; ++i) {
    ground += power[i];
  }
  ground /= 60.0;
  CHECK(10.0 * std::log10(power[static_cast<std::size_t>(peak)] / ground) > 15.0);
}

TEST_CASE("Backprojection kernels agree bit for bit", "[coherent][backprojection]") {
  const RadarParams radar = theory_radar();
  const LinearTrajectory track({-4000.0, 0.0, 3000.0}, {0.0, 100.0, 0.0}, -10.0, 10.0);
  const std::vector<Scatterer> s{{{0.0, 0.0, 0.0}, 1.0}, {{3.0, 2.0, 0.0}, 0.5}};
  const EchoWindow w{-0.7, 0.7, 4950.0, 5050.0};
  const CompressedData comp = range_compress(synthesize_raw(s, track, radar, w), 4);
  PixelGrid grid;
  grid.width = 40;
  grid.height = 30;
  grid.origin = {-5.0, -5.0, 0.0};
  grid.step_x = {0.25, 0.0, 0.0};
  grid.step_y = {0.0, 0.25, 0.0};
  BackprojectionConfig naive;
  naive.kernel = BackprojectionKernel::Naive;
  BackprojectionConfig blocked;
  blocked.pixel_block = 37;
  blocked.pulse_block = 5;
  const auto a = backproject(comp, grid, naive);
  const auto b = backproject(comp, grid, blocked);
  REQUIRE(a.data.size() == b.data.size());
  for (std::size_t i = 0; i < a.data.size(); ++i) {
    REQUIRE(a.data[i] == b.data[i]);
  }
  // The brightest pixel is at the stronger scatterer (origin = pixel 20, 20).
  const auto peak = std::max_element(a.data.begin(), a.data.end(),
                                     [](auto x, auto y) { return std::abs(x) < std::abs(y); }) -
                    a.data.begin();
  CHECK(static_cast<std::size_t>(peak) == 20 * grid.width + 20);
}
