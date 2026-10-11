#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <limits>
#include <vector>

#include "backscatter/image/geocode.hpp"
#include "backscatter/math/constants.hpp"
#include "backscatter/scene/heightmap.hpp"
#include "backscatter/scene/scene.hpp"
#include "backscatter/sensor/passes.hpp"
#include "backscatter/sensor/sar_geometry.hpp"
#include "backscatter/trace/geometric.hpp"

using namespace bsar;
using Catch::Approx;

TEST_CASE("format_utc inverts parse_utc", "[sensor][passes]") {
  for (const char* s : {"2000-01-01T00:00:00.000000", "2024-02-29T23:59:58.250000",
                        "2026-10-10T12:34:56.789000", "1999-12-31T06:00:00.000000"}) {
    CHECK(format_utc(parse_utc(s)) == s);
  }
}

TEST_CASE("Pass finder on a synthetic polar orbit", "[sensor][passes]") {
  // Circular polar orbit (Earth rotation ignored) whose plane lies 4.5 deg of
  // longitude west of the target: moving north it sees the target on its right.
  const Geodetic target{46.5, 7.9, 1000.0};
  const double radius = 7.071e6;
  const double omega = std::sqrt(3.986004418e14 / (radius * radius * radius));
  const double lon = deg_to_rad(target.lon_deg - 4.5);
  std::vector<StateVector> svs;
  for (int k = 0; k <= 600; ++k) {  // 100 minutes, one revolution
    const double t = 10.0 * k;
    const double a = omega * t;
    svs.push_back({t,
                   {radius * std::cos(a) * std::cos(lon), radius * std::cos(a) * std::sin(lon),
                    radius * std::sin(a)},
                   {-radius * omega * std::sin(a) * std::cos(lon),
                    -radius * omega * std::sin(a) * std::sin(lon), radius * omega * std::cos(a)}});
  }
  const auto passes = find_passes(Orbit(svs), target);
  REQUIRE(passes.size() == 1);
  const OrbitPass& p = passes.front();
  CHECK(p.ascending);
  CHECK(p.side == LookSide::Right);
  CHECK(p.incidence_deg > 20.0);
  CHECK(p.incidence_deg < 40.0);
  CHECK(p.range > 7.0e5);
  CHECK(p.range < 9.0e5);
  // Closest approach is near the time the satellite crosses latitude 46.5.
  CHECK(p.t == Approx(deg_to_rad(46.3) / omega).margin(60.0));
}

TEST_CASE("Orbits keep precise time when the epoch is far from 2000", "[sensor][passes]") {
  // Same geometry as above, but stamped in August 2024 (t ~ 7.8e8 s), where an
  // absolute double only resolves ~1.2e-7 s. Orbit time keeps sub-ns precision
  // and the zero-Doppler solve converges.
  const Geodetic target{46.5, 7.9, 1000.0};
  const double epoch = parse_utc("2024-08-15T05:00:00");
  const double radius = 7.071e6;
  const double omega = std::sqrt(3.986004418e14 / (radius * radius * radius));
  const double lon = deg_to_rad(target.lon_deg - 4.5);
  std::vector<StateVector> svs;
  for (int k = 0; k <= 600; ++k) {
    const double t = 10.0 * k;
    const double a = omega * t;
    svs.push_back({epoch + t,
                   {radius * std::cos(a) * std::cos(lon), radius * std::cos(a) * std::sin(lon),
                    radius * std::sin(a)},
                   {-radius * omega * std::sin(a) * std::cos(lon),
                    -radius * omega * std::sin(a) * std::sin(lon), radius * omega * std::cos(a)}});
  }
  const Orbit orbit(svs);
  CHECK(orbit.epoch() == epoch);
  CHECK(orbit.t_begin() == 0.0);
  const auto passes = find_passes(orbit, target);
  REQUIRE(passes.size() == 1);
  CHECK(passes.front().t - epoch == Approx(deg_to_rad(46.3) / omega).margin(60.0));
  const OrbitTrajectory traj(orbit, EnuFrame(target));
  const ZeroDopplerSolution zd =
      solve_zero_doppler(traj, Vec3d{}, passes.front().t - orbit.epoch() + 3.0);
  CHECK(zd.converged);
  CHECK(std::abs(dot(traj.velocity(zd.t), traj.position(zd.t))) < 1e-3 * 7.5e3);
  // A window keeps absolute alignment.
  const Orbit w = orbit.window(epoch + 3000.0, 600.0);
  CHECK(w.epoch() == epoch + 2400.0);
  CHECK(length(w.position(600.0) - orbit.position(3000.0)) < 1e-6);
}

TEST_CASE("Geocoded masks match the closed-form ridge geometry", "[geocode]") {
  // Ridge along y on a 10 m grid: west face 50 deg (> 35 deg incidence, so in
  // layover), east face 70 deg (> 55 deg, so it and the ground behind it are
  // in shadow out to H tan(35) = 210 m east of the crest).
  constexpr double kH = 300.0;
  const double west_run = kH / std::tan(deg_to_rad(50.0));
  const double east_run = kH / std::tan(deg_to_rad(70.0));
  Heightmap hm;
  hm.nx = 121;
  hm.ny = 61;
  hm.dx = hm.dy = 10.0;
  hm.x0 = -600.0;
  hm.y0 = -300.0;
  hm.z.resize(hm.nx * hm.ny);
  for (std::size_t j = 0; j < hm.ny; ++j) {
    for (std::size_t i = 0; i < hm.nx; ++i) {
      const double x = hm.x0 + static_cast<double>(i) * hm.dx;
      double z = 0.0;
      if (x > -west_run && x <= 0.0) {
        z = kH * (1.0 + x / west_run);
      } else if (x > 0.0 && x < east_run) {
        z = kH * (1.0 - x / east_run);
      }
      hm.at(i, j) = static_cast<float>(z);
    }
  }
  Scene scene;
  scene.mesh = heightmap_to_mesh(hm, scene.material_id("dry_soil"));
  scene.build();
  const auto track = make_side_looking_track(scene.bounds(), 693e3, 7600.0, 35.0, LookSide::Right);
  GeometricConfig cfg;
  cfg.range_spacing = 2.0;
  cfg.azimuth_spacing = 10.0;
  cfg.max_bounces = 1;
  const GeometricImage img = render_geometric(scene, track, cfg);

  std::vector<Vec3d> points(hm.nx * hm.ny);
  for (std::size_t j = 0; j < hm.ny; ++j) {
    for (std::size_t i = 0; i < hm.nx; ++i) {
      points[j * hm.nx + i] = {hm.x0 + static_cast<double>(i) * hm.dx,
                               hm.y0 + static_cast<double>(j) * hm.dy, hm.at(i, j)};
    }
  }
  points[0].z = std::numeric_limits<double>::quiet_NaN();  // a no-data sample
  const GeocodedProducts geo = geocode(img, track, scene, points, hm.nx, hm.ny);

  const std::size_t row = hm.ny / 2;
  auto col = [&](double x) { return static_cast<std::size_t>(std::lround((x - hm.x0) / hm.dx)); };
  CHECK(geo.valid.at(0, 0) == 0);
  CHECK(geo.valid.at(col(-400.0), row) == 1);
  // Middle of the west face: layover, lit.
  CHECK(geo.layover.at(col(-130.0), row) == 1);
  CHECK(geo.shadow.at(col(-130.0), row) == 0);
  // Passive layover: flat ground in front of the face whose range bins it shares.
  CHECK(geo.layover.at(col(-400.0), row) == 1);
  // Ground behind the ridge, between the east foot (109 m) and 210 m: shadow.
  CHECK(geo.shadow.at(col(160.0), row) == 1);
  // Beyond the shadow, and in front of the passive-layover zone (the crest has
  // the same slant range as ground H / tan(35) = 428 m in front of it): neither.
  CHECK(geo.shadow.at(col(300.0), row) == 0);
  CHECK(geo.layover.at(col(-550.0), row) == 0);
  CHECK(geo.shadow.at(col(-550.0), row) == 0);
  // Flat ground there renders beta0 = sigma0 / sin(theta).
  const Vec3d p = track.position(0.0);
  const double inc = std::atan2(-550.0 - p.x, p.z);
  const double sigma0 =
      FacetScatteringModel{}.backscatter(materials::dry_soil(), std::cos(inc), Polarization::VV);
  CHECK(geo.beta0.at(col(-550.0), row) == Approx(sigma0 / std::sin(inc)).epsilon(0.1));
}
