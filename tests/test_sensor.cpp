#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <string>

#include "backscatter/math/constants.hpp"
#include "backscatter/sensor/geodesy.hpp"
#include "backscatter/sensor/orbit.hpp"
#include "backscatter/sensor/radar.hpp"
#include "backscatter/sensor/sar_geometry.hpp"
#include "backscatter/sensor/trajectory.hpp"

using namespace bsar;
using Catch::Approx;

TEST_CASE("Geodetic <-> ECEF round trip", "[sensor][geodesy]") {
  const Vec3d equator = geodetic_to_ecef({0.0, 0.0, 0.0});
  CHECK(equator.x == Approx(wgs84::kSemiMajor));
  CHECK(std::abs(equator.y) < 1e-6);
  const Vec3d pole = geodetic_to_ecef({90.0, 0.0, 0.0});
  CHECK(pole.z == Approx(wgs84::kSemiMinor));

  for (const Geodetic g : {Geodetic{46.5, 7.9, 3500.0}, Geodetic{-33.9, 151.2, 20.0},
                           Geodetic{89.9, -45.0, 100.0}, Geodetic{0.1, 179.9, -50.0}}) {
    const Geodetic back = ecef_to_geodetic(geodetic_to_ecef(g));
    CHECK(back.lat_deg == Approx(g.lat_deg).margin(1e-9));
    CHECK(back.lon_deg == Approx(g.lon_deg).margin(1e-9));
    CHECK(back.height == Approx(g.height).margin(1e-4));
  }
}

TEST_CASE("ENU frame axes", "[sensor][geodesy]") {
  const EnuFrame frame({46.5, 7.9, 500.0});
  const Vec3d origin = frame.geodetic_to_enu({46.5, 7.9, 500.0});
  CHECK(length(origin) < 1e-6);
  const Vec3d up = frame.geodetic_to_enu({46.5, 7.9, 1500.0});
  CHECK(up.z == Approx(1000.0).margin(1e-6));
  CHECK(std::abs(up.x) < 1e-6);
  const Vec3d north = frame.geodetic_to_enu({46.51, 7.9, 500.0});
  CHECK(north.y > 1000.0);
  CHECK(std::abs(north.x) < 1e-6);
  const Vec3d east = frame.geodetic_to_enu({46.5, 7.91, 500.0});
  CHECK(east.x > 700.0);
  const Vec3d v{1.0, 2.0, 3.0};
  const Vec3d back = frame.rotate_to_ecef(frame.rotate_to_enu(v));
  CHECK(length(back - v) < 1e-12);
}

TEST_CASE("UTC parsing", "[sensor][orbit]") {
  CHECK(parse_utc("2000-01-01T00:00:00") == 0.0);
  CHECK(parse_utc("UTC=2000-01-02T00:00:01.5") == Approx(86401.5));
  // 2024 is a leap year: 2024-03-01 is 31 + 29 days after 2024-01-01.
  CHECK(parse_utc("2024-03-01T00:00:00") - parse_utc("2024-01-01T00:00:00") == 60.0 * 86400.0);
  CHECK_THROWS(parse_utc("yesterday"));
}

TEST_CASE("Hermite orbit interpolation on a circular orbit", "[sensor][orbit]") {
  // Circular orbit at Sentinel-1 radius, state vectors every 10 s.
  const double radius = 7.071e6;
  const double omega = std::sqrt(3.986004418e14 / (radius * radius * radius));
  auto pos = [&](double t) {
    return Vec3d{radius * std::cos(omega * t), radius * std::sin(omega * t) * 0.6,
                 radius * std::sin(omega * t) * 0.8};
  };
  auto vel = [&](double t) {
    return Vec3d{-radius * omega * std::sin(omega * t), radius * omega * std::cos(omega * t) * 0.6,
                 radius * omega * std::cos(omega * t) * 0.8};
  };
  std::vector<StateVector> svs;
  for (int k = 0; k < 30; ++k) {
    const double t = 10.0 * k;
    svs.push_back({t, pos(t), vel(t)});
  }
  const Orbit orbit(svs);
  double max_pos_err = 0.0;
  double max_vel_err = 0.0;
  for (double t = 0.0; t <= 290.0; t += 0.37) {
    max_pos_err = std::max(max_pos_err, length(orbit.position(t) - pos(t)));
    max_vel_err = std::max(max_vel_err, length(orbit.velocity(t) - vel(t)));
  }
  CHECK(max_pos_err < 1e-3);  // < 1 mm
  CHECK(max_vel_err < 1e-2);  // < 1 cm/s
}

TEST_CASE("EOF orbit file parsing", "[sensor][orbit]") {
  const std::string xml = R"(<?xml version="1.0"?>
<Earth_Explorer_File><Data_Block type="xml"><List_of_OSVs count="2">
<OSV>
  <TAI>TAI=2024-01-01T23:00:37.000000</TAI>
  <UTC>UTC=2024-01-01T23:00:00.000000</UTC>
  <UT1>UT1=2024-01-01T23:00:00.012345</UT1>
  <Absolute_Orbit>+51893</Absolute_Orbit>
  <X unit="m">-1234567.891</X>
  <Y unit="m">2345678.912</Y>
  <Z unit="m">6543210.123</Z>
  <VX unit="m/s">1000.5</VX>
  <VY unit="m/s">-7000.25</VY>
  <VZ unit="m/s">2500.125</VZ>
  <Quality>NOMINAL</Quality>
</OSV>
<OSV>
  <TAI>TAI=2024-01-01T23:00:47.000000</TAI>
  <UTC>UTC=2024-01-01T23:00:10.000000</UTC>
  <UT1>UT1=2024-01-01T23:00:10.012345</UT1>
  <Absolute_Orbit>+51893</Absolute_Orbit>
  <X unit="m">-1224562.886</X>
  <Y unit="m">2275676.410</Y>
  <Z unit="m">6568211.373</Z>
  <VX unit="m/s">1000.5</VX>
  <VY unit="m/s">-7000.25</VY>
  <VZ unit="m/s">2500.125</VZ>
  <Quality>NOMINAL</Quality>
</OSV>
</List_of_OSVs></Data_Block></Earth_Explorer_File>)";
  const Orbit orbit = parse_eof(xml);
  REQUIRE(orbit.state_vectors().size() == 2);
  const StateVector& sv = orbit.state_vectors()[0];
  CHECK(sv.t == Approx(parse_utc("2024-01-01T23:00:00")));
  CHECK(sv.position.x == Approx(-1234567.891));
  CHECK(sv.velocity.z == Approx(2500.125));
  CHECK(orbit.t_end() - orbit.t_begin() == Approx(10.0));
  // Endpoints are reproduced exactly.
  CHECK(length(orbit.position(orbit.t_begin()) - sv.position) < 1e-6);
}

TEST_CASE("Zero-Doppler solver", "[sensor][geometry]") {
  const LinearTrajectory track({0.0, -5000.0, 7000.0}, {0.0, 120.0, 0.0}, -1000.0, 1000.0);
  const Vec3d target{4000.0, 1234.0, 50.0};
  const auto sol = solve_zero_doppler(track, target, -50.0);
  REQUIRE(sol.converged);
  CHECK(sol.t == Approx((1234.0 + 5000.0) / 120.0).epsilon(1e-12));
  CHECK(sol.range == Approx(std::hypot(4000.0, 7000.0 - 50.0)));

  // Curved track: verify the defining property v . (p - x) = 0.
  const double radius = 7.0e6;
  const double omega = 1e-3;
  struct Circle final : Trajectory {
    double r, w;
    Circle(double r_, double w_) : r(r_), w(w_) {}
    Vec3d position(double t) const override {
      return {r * std::cos(w * t), r * std::sin(w * t), 0};
    }
    Vec3d velocity(double t) const override {
      return {-r * w * std::sin(w * t), r * w * std::cos(w * t), 0};
    }
    double t_begin() const override { return -1000; }
    double t_end() const override { return 1000; }
  } circle(radius, omega);
  const Vec3d ground{6.3e6, 3.0e5, 1.0e5};
  const auto s2 = solve_zero_doppler(circle, ground, 0.0);
  REQUIRE(s2.converged);
  const double doppler = dot(circle.velocity(s2.t), circle.position(s2.t) - ground);
  CHECK(std::abs(doppler) < 1e-3 * length(circle.velocity(s2.t)));
}

TEST_CASE("Side-looking track geometry", "[sensor][geometry]") {
  Aabb box;
  box.expand(Vec3f{-100, -100, 0});
  box.expand(Vec3f{100, 100, 0});
  const auto right = make_side_looking_track(box, 1000.0, 50.0, 30.0, LookSide::Right);
  const Vec3d p = right.position(0.0);
  CHECK(p.x < 0.0);  // west of the scene, looking east
  CHECK(std::atan2(-p.x, p.z) == Approx(deg_to_rad(30.0)));
  const auto left = make_side_looking_track(box, 1000.0, 50.0, 30.0, LookSide::Left);
  CHECK(left.position(0.0).x > 0.0);
}

TEST_CASE("Chirp sweeps the requested band", "[sensor][radar]") {
  RadarParams r;
  r.bandwidth = 20e6;
  r.pulse_duration = 10e-6;
  r.sample_rate = 40e6;
  const auto chirp = make_chirp(r);
  CHECK(chirp.size() == 400);
  // Instantaneous frequency from the phase difference at both ends.
  const double f_start = std::arg(chirp[1] * std::conj(chirp[0])) * r.sample_rate / (2 * kPi);
  const double f_end = std::arg(chirp[399] * std::conj(chirp[398])) * r.sample_rate / (2 * kPi);
  CHECK(f_start == Approx(-10e6).epsilon(0.01));
  CHECK(f_end == Approx(10e6).epsilon(0.01));
  CHECK(slant_range_resolution(150e6) == Approx(0.99930819333));
  RadarParams bad = r;
  bad.sample_rate = 10e6;
  CHECK_THROWS(make_chirp(bad));
}
