#include <algorithm>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>

#include "backscatter/eval/theory_checks.hpp"
#include "backscatter/math/constants.hpp"
#include "backscatter/scene/scene.hpp"
#include "backscatter/sensor/sar_geometry.hpp"
#include "backscatter/trace/geometric.hpp"
#include "backscatter/trace/hillshade.hpp"

using namespace bsar;
using Catch::Approx;

TEST_CASE("Flat plane renders beta0 = sigma0 / sin(incidence)", "[geometric]") {
  const Scene scene = make_plane_scene(2000.0);
  constexpr double kIncidence = 40.0;
  const auto track =
      make_side_looking_track(scene.bounds(), 693e3, 7600.0, kIncidence, LookSide::Right);
  GeometricConfig cfg;
  cfg.range_spacing = 10.0;
  cfg.azimuth_spacing = 50.0;
  cfg.max_bounces = 2;
  const GeometricImage img = render_geometric(scene, track, cfg);

  REQUIRE(img.intensity.width > 50);
  const std::size_t line = img.intensity.height / 2;
  const std::size_t mid = img.intensity.width / 2;
  const double measured = img.intensity.at(mid, line);
  // Incidence at the centre bin's ground position.
  const Vec3d p = track.position(img.t_start + static_cast<double>(line) * img.line_interval);
  const double r = img.near_range + (static_cast<double>(mid) + 0.5) * img.range_spacing;
  const double ground = std::sqrt(r * r - p.z * p.z);
  const double inc = std::atan2(ground, p.z);
  const FacetScatteringModel model;
  const double sigma0 = model.backscatter(materials::dry_soil(), std::cos(inc), Polarization::VV);
  CHECK(measured == Approx(sigma0 / std::sin(inc)).epsilon(0.05));
  // A flat plane has neither layover nor shadow, and no multi-bounce.
  CHECK(std::count(img.layover.data.begin(), img.layover.data.end(), 1) == 0);
  CHECK(std::count(img.shadow.data.begin(), img.shadow.data.end(), 1) == 0);
  CHECK(*std::max_element(img.bounce[1].data.begin(), img.bounce[1].data.end()) == 0.0f);
}

TEST_CASE("Layover and shadow match the closed-form slope conditions", "[geometric][theory]") {
  for (const LayoverShadowCase& c : run_layover_shadow_checks()) {
    INFO("west " << c.west_slope_deg << " east " << c.east_slope_deg << " layover "
                 << c.layover_fraction << " shadow " << c.shadow_fraction);
    CHECK(c.pass());
  }
}

TEST_CASE("Dihedral double bounce lands at the corner range", "[geometric][theory]") {
  const DihedralReport rep = run_dihedral_check();
  INFO("predicted " << rep.predicted_range << " measured " << rep.measured_range);
  CHECK(std::abs(rep.measured_range - rep.predicted_range) <= rep.range_spacing);
  CHECK(rep.double_to_single_db > 3.0);  // the corner line outshines single bounce
}

TEST_CASE("Rendering is bit-identical across thread counts", "[geometric][coherent][theory]") {
  const DeterminismReport rep = run_determinism_check(8);
  CHECK(rep.geometric_identical);
  CHECK(rep.coherent_identical);
}

TEST_CASE("Hillshade of a tilted plane matches the analytic value", "[hillshade]") {
  // z = 0.3 x: normal = (-0.3, 0, 1) / |.|; sun from the west at 45 deg.
  Scene scene;
  const auto m = scene.material_id("dry_soil");
  const auto a = scene.mesh.add_vertex({-100, -100, -30});
  const auto b = scene.mesh.add_vertex({100, -100, 30});
  const auto c = scene.mesh.add_vertex({100, 100, 30});
  const auto d = scene.mesh.add_vertex({-100, 100, -30});
  scene.mesh.add_triangle(a, b, c, m);
  scene.mesh.add_triangle(a, c, d, m);
  scene.build();
  HillshadeConfig cfg;
  cfg.sun_azimuth_deg = 270.0;
  cfg.sun_elevation_deg = 45.0;
  cfg.pixel_size = 5.0;
  const Image<float> img = render_hillshade(scene, cfg);
  const Vec3d n = normalize(Vec3d{-0.3, 0.0, 1.0});
  const Vec3d sun{-std::cos(deg_to_rad(45.0)), 0.0, std::sin(deg_to_rad(45.0))};
  const double expected = dot(n, sun);
  CHECK(img.at(img.width / 2, img.height / 2) == Approx(expected).epsilon(1e-5));
  CHECK(img.at(1, 1) == Approx(expected).epsilon(1e-5));
}

TEST_CASE("Sensor beyond float range still bins correctly", "[geometric]") {
  // A point-like scene seen from 900 km: the binned range must match the
  // double-precision distance to within a bin.
  const Scene scene = make_dihedral_scene(60.0, 4.0, 4.0, 4.0);
  const auto track = make_side_looking_track(scene.bounds(), 800e3, 7500.0, 30.0, LookSide::Left);
  GeometricConfig cfg;
  cfg.range_spacing = 0.5;
  cfg.azimuth_spacing = 2.0;
  cfg.max_bounces = 1;
  const GeometricImage img = render_geometric(scene, track, cfg);
  const std::size_t line = img.intensity.height / 2;
  const float* row = img.intensity.row(line);
  std::size_t first = 0;
  while (first < img.intensity.width && row[first] == 0.0f) {
    ++first;
  }
  REQUIRE(first < img.intensity.width);
  // The nearest surface to a left-looking sensor (east of the scene) is the
  // east edge of the ground (it is ~10 m closer along the line of sight than
  // the roof's east edge).
  const Vec3d p = track.position(img.t_start + static_cast<double>(line) * img.line_interval);
  const double nearest = length(p - Vec3d{30.0, p.y, 0.0});
  const double first_range = img.near_range + static_cast<double>(first) * img.range_spacing;
  CHECK(std::abs(first_range - nearest) < 2.0 * cfg.range_spacing);
}
