#include "backscatter/eval/theory_checks.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>

#include "backscatter/image/backprojection.hpp"
#include "backscatter/scene/scene.hpp"
#include "backscatter/sensor/sar_geometry.hpp"
#include "backscatter/trace/coherent.hpp"
#include "backscatter/trace/geometric.hpp"
#include "backscatter/util/parallel.hpp"

namespace bsar {
namespace {

constexpr double kOrbitAltitude = 693e3;  // Sentinel-1 nominal
constexpr double kOrbitSpeed = 7600.0;

/// Straight track heading north at x = -ground_range, passing y = 0 at t = 0.
LinearTrajectory airborne_track(const CoherentGeometry& g) {
  return LinearTrajectory({-g.ground_range, 0.0, g.altitude}, {0.0, g.speed, 0.0}, -1e4, 1e4);
}

std::vector<float> amplitude_cut(const CompressedData& data, const Vec3d& center,
                                 const Vec3d& direction, double spacing, std::size_t n,
                                 unsigned threads) {
  PixelGrid grid;
  grid.width = n;
  grid.height = 1;
  const std::size_t center_index = n / 2;  // n is odd: the centre sample
  grid.origin = center - direction * (spacing * static_cast<double>(center_index));
  grid.step_x = direction * spacing;
  BackprojectionConfig cfg;
  cfg.threads = threads;
  const auto img = backproject(data, grid, cfg);
  std::vector<float> amp(n);
  for (std::size_t i = 0; i < n; ++i) {
    amp[i] = std::abs(img.data[i]);
  }
  return amp;
}

template <typename T>
bool same_bits(const std::vector<T>& a, const std::vector<T>& b) {
  return a.size() == b.size() && std::memcmp(a.data(), b.data(), a.size() * sizeof(T)) == 0;
}

}  // namespace

RadarParams theory_radar() {
  RadarParams r;
  r.carrier_frequency = 5.405e9;
  r.bandwidth = 100e6;
  r.pulse_duration = 5e-6;
  r.sample_rate = 150e6;
  r.antenna_length = 2.0;
  r.prf = 200.0;  // 0.5 m pulse spacing at 100 m/s = L / 4
  r.beam = BeamPattern::Uniform;
  return r;
}

PointTargetReport run_point_target_check(const RadarParams& radar, const CoherentGeometry& geom,
                                         unsigned threads) {
  const LinearTrajectory track = airborne_track(geom);
  const Vec3d target{0.0, 0.0, 0.0};
  const double r0 = length(target - track.position(0.0));
  const double half_angle =
      (radar.beam == BeamPattern::Uniform ? 0.5 : 1.0) * radar.wavelength() / radar.antenna_length;
  const double half_aperture = r0 * std::tan(half_angle);
  const double res_r = slant_range_resolution(radar.bandwidth);
  const double res_a = stripmap_azimuth_resolution(radar.antenna_length);

  EchoWindow w;
  w.t_start = -half_aperture / geom.speed - 2.0 / radar.prf;
  w.t_end = half_aperture / geom.speed + 2.0 / radar.prf;
  w.near_range = r0 - 40.0 * res_r;
  w.far_range = std::hypot(r0, half_aperture) + 40.0 * res_r;

  const std::vector<Scatterer> targets{{target, 1.0}};
  const RawData raw = synthesize_raw(targets, track, radar, w, threads);
  const CompressedData comp = range_compress(raw, 8, threads);

  constexpr std::size_t kCut = 961;  // +/- 24 resolution cells at 1/20 cell spacing
  const Vec3d los = normalize(target - track.position(0.0));
  const auto range_amp = amplitude_cut(comp, target, los, res_r / 20.0, kCut, threads);
  const auto az_amp = amplitude_cut(comp, target, {0.0, 1.0, 0.0}, res_a / 20.0, kCut, threads);

  PointTargetReport rep;
  rep.pulses = raw.num_pulses;
  rep.range = analyze_irf(range_amp, res_r / 20.0);
  rep.azimuth = analyze_irf(az_amp, res_a / 20.0);
  rep.expected_range_3db = kSinc3dbFactor * res_r;
  rep.expected_azimuth_3db = kSinc3dbFactor * res_a;
  rep.range_error_pct = 100.0 * (rep.range.resolution_3db / rep.expected_range_3db - 1.0);
  rep.azimuth_error_pct = 100.0 * (rep.azimuth.resolution_3db / rep.expected_azimuth_3db - 1.0);
  return rep;
}

SpeckleReport run_speckle_check(const RadarParams& radar, const CoherentGeometry& geom,
                                double patch_size, double density, std::uint64_t seed,
                                unsigned threads) {
  const Scene scene = make_plane_scene(patch_size);
  const LinearTrajectory track = airborne_track(geom);
  ScattererConfig sc;
  sc.density = density;
  sc.seed = seed;
  const FacetScatteringModel model;
  const auto scatterers = sample_scatterers(scene, track.position(0.0), sc, model);
  const EchoWindow w = auto_echo_window(scene.bounds(), track, radar);
  const RawData raw = synthesize_raw(scatterers, track, radar, w, threads);
  const CompressedData comp = range_compress(raw, 8, threads);

  // Two resolution cells between samples in both directions.
  const double incidence = std::atan2(geom.ground_range, geom.altitude);
  const double step_ground = 2.0 * slant_range_resolution(radar.bandwidth) / std::sin(incidence);
  const double step_az = 2.0 * stripmap_azimuth_resolution(radar.antenna_length);
  const double interior = 0.7 * patch_size;  // stay clear of the patch edges
  PixelGrid grid;
  grid.width = static_cast<std::size_t>(interior / step_ground) + 1;
  grid.height = static_cast<std::size_t>(interior / step_az) + 1;
  grid.origin = {-0.5 * interior, -0.5 * interior, 0.0};
  grid.step_x = {step_ground, 0.0, 0.0};
  grid.step_y = {0.0, step_az, 0.0};
  BackprojectionConfig cfg;
  cfg.threads = threads;
  const auto img = backproject(comp, grid, cfg);

  std::vector<double> intensity(img.size());
  for (std::size_t i = 0; i < img.size(); ++i) {
    intensity[i] = std::norm(std::complex<double>(img.data[i]));
  }
  SpeckleReport rep;
  rep.stats = speckle_statistics(intensity);
  rep.scatterers = scatterers.size();
  rep.pulses = raw.num_pulses;
  return rep;
}

std::vector<LayoverShadowCase> run_layover_shadow_checks(unsigned threads) {
  constexpr double kIncidence = 35.0;
  const std::vector<std::pair<double, double>> slopes{{45.0, 20.0}, {25.0, 20.0}, {25.0, 65.0},
                                                      {45.0, 65.0}, {30.0, 50.0}, {40.0, 60.0}};
  std::vector<LayoverShadowCase> out;
  for (const auto& [west, east] : slopes) {
    const Scene scene = make_ridge_scene(west, east, 300.0, 3000.0);
    const auto track = make_side_looking_track(scene.bounds(), kOrbitAltitude, kOrbitSpeed,
                                               kIncidence, LookSide::Right);
    GeometricConfig cfg;
    cfg.range_spacing = 5.0;
    cfg.azimuth_spacing = 100.0;
    cfg.max_bounces = 1;
    cfg.threads = threads;
    const GeometricImage img = render_geometric(scene, track, cfg);

    LayoverShadowCase c;
    c.incidence_deg = kIncidence;
    c.west_slope_deg = west;
    c.east_slope_deg = east;
    c.expect_layover = west > kIncidence;
    c.expect_shadow = east > 90.0 - kIncidence;
    const auto n = static_cast<double>(img.layover.size());
    c.layover_fraction =
        static_cast<double>(std::count(img.layover.data.begin(), img.layover.data.end(), 1)) / n;
    c.shadow_fraction =
        static_cast<double>(std::count(img.shadow.data.begin(), img.shadow.data.end(), 1)) / n;
    out.push_back(c);
  }
  return out;
}

DihedralReport run_dihedral_check(unsigned threads) {
  constexpr double kWidth = 40.0;
  const Scene scene = make_dihedral_scene(400.0, kWidth, kWidth, 30.0);
  const auto track =
      make_side_looking_track(scene.bounds(), kOrbitAltitude, kOrbitSpeed, 35.0, LookSide::Right);
  GeometricConfig cfg;
  cfg.range_spacing = 1.0;
  cfg.azimuth_spacing = 4.0;
  cfg.max_bounces = 3;
  cfg.threads = threads;
  const GeometricImage img = render_geometric(scene, track, cfg);

  // Azimuth line through the middle of the sensor-facing (west) wall.
  const auto line = static_cast<std::size_t>(std::lround(-img.t_start / img.line_interval));
  const double t = img.t_start + static_cast<double>(line) * img.line_interval;
  const Vec3d p = track.position(t);
  const Vec3d corner{-0.5 * kWidth, p.y, 0.0};

  const float* dbl = img.bounce[1].row(line);
  const float* sgl = img.bounce[0].row(line);
  const auto peak2 = std::max_element(dbl, dbl + img.bounce[1].width) - dbl;
  const float peak1 = *std::max_element(sgl, sgl + img.bounce[0].width);

  DihedralReport rep;
  rep.predicted_range = length(p - corner);
  rep.measured_range = img.near_range + (static_cast<double>(peak2) + 0.5) * img.range_spacing;
  rep.range_spacing = img.range_spacing;
  rep.double_to_single_db = 10.0 * std::log10(dbl[peak2] / peak1);
  return rep;
}

DeterminismReport run_determinism_check(unsigned threads) {
  DeterminismReport rep;
  rep.threads_compared = std::max(2u, resolve_threads(threads));

  {
    const Scene scene = make_synthetic_scene("city");
    const auto track =
        make_side_looking_track(scene.bounds(), kOrbitAltitude, kOrbitSpeed, 35.0, LookSide::Right);
    GeometricConfig cfg;
    cfg.range_spacing = 2.0;
    cfg.azimuth_spacing = 8.0;
    cfg.threads = 1;
    const GeometricImage a = render_geometric(scene, track, cfg);
    cfg.threads = rep.threads_compared;
    const GeometricImage b = render_geometric(scene, track, cfg);
    bool same = same_bits(a.intensity.data, b.intensity.data) &&
                same_bits(a.layover.data, b.layover.data) &&
                same_bits(a.shadow.data, b.shadow.data);
    for (std::size_t k = 0; k < a.bounce.size(); ++k) {
      same = same && same_bits(a.bounce[k].data, b.bounce[k].data);
    }
    rep.geometric_identical = same;
  }
  {
    const RadarParams radar = theory_radar();
    const CoherentGeometry geom;
    const Scene scene = make_plane_scene(30.0);
    const LinearTrajectory track = airborne_track(geom);
    ScattererConfig sc;
    sc.density = 2.0;
    auto scatterers = sample_scatterers(scene, track.position(0.0), sc, FacetScatteringModel{});
    scatterers.push_back({{0.0, 0.0, 0.0}, 5.0});
    const EchoWindow w = auto_echo_window(scene.bounds(), track, radar);
    auto run = [&](unsigned th) {
      const RawData raw = synthesize_raw(scatterers, track, radar, w, th);
      const CompressedData comp = range_compress(raw, 4, th);
      PixelGrid grid;
      grid.width = grid.height = 48;
      grid.origin = {-12.0, -12.0, 0.0};
      grid.step_x = {0.5, 0.0, 0.0};
      grid.step_y = {0.0, 0.5, 0.0};
      BackprojectionConfig cfg;
      cfg.threads = th;
      return backproject(comp, grid, cfg).data;
    };
    rep.coherent_identical = same_bits(run(1), run(rep.threads_compared));
  }
  return rep;
}

}  // namespace bsar
