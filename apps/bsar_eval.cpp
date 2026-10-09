// bsar_eval: validation checks and image comparison metrics.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>

#include "backscatter/eval/metrics.hpp"
#include "backscatter/eval/theory_checks.hpp"
#include "backscatter/image/io.hpp"
#include "backscatter/sensor/sar_geometry.hpp"

#include "cli.hpp"

using namespace bsar;

namespace {

constexpr const char* kUsage = R"(usage: bsar_eval <command> [options]

  theory [--threads N]        run the self-contained validation suite (point
                              target, speckle, layover/shadow, dihedral,
                              determinism); prints a Markdown report and exits
                              non-zero if any check fails
  irf <slc.npy> --dx M --dy M impulse response metrics through the brightest
                              pixel (column spacing dx, row spacing dy)
  speckle <image.npy>         ENL and exponential goodness of fit; complex or
                              amplitude input is squared to intensity unless
                              --intensity is given
  iou <a.npy> <b.npy>         intersection over union of two masks
)";

const char* mark(bool ok) { return ok ? "pass" : "**FAIL**"; }

int run_theory(const cli::Args& args) {
  const auto threads = static_cast<unsigned>(args.num("threads", 0));
  args.reject_unknown();
  bool all = true;
  const RadarParams radar = theory_radar();
  const CoherentGeometry geom;

  std::printf("## Theory checks\n\n");
  std::printf(
      "System: %.3f GHz carrier, %.0f MHz chirp (%.1f us), %.0f MHz sampling, %.1f m "
      "antenna, unweighted; airborne track at %.0f m altitude, %.0f m ground range.\n\n",
      radar.carrier_frequency * 1e-9, radar.bandwidth * 1e-6, radar.pulse_duration * 1e6,
      radar.sample_rate * 1e-6, radar.antenna_length, geom.altitude, geom.ground_range);

  const PointTargetReport pt = run_point_target_check(radar, geom, threads);
  const bool range_ok = std::abs(pt.range_error_pct) < 5.0;
  const bool az_ok = std::abs(pt.azimuth_error_pct) < 5.0;
  const bool pslr_ok = std::abs(pt.range.pslr_db - pt.expected_pslr_db) < 1.0 &&
                       std::abs(pt.azimuth.pslr_db - pt.expected_pslr_db) < 1.0;
  all = all && range_ok && az_ok && pslr_ok;
  std::printf("### Point target (%zu pulses)\n\n", pt.pulses);
  std::printf("| Metric | Measured | Theory | Error | Criterion | Result |\n");
  std::printf("|---|---|---|---|---|---|\n");
  std::printf("| Range -3 dB width | %.4f m | %.4f m | %+.2f%% | within 5%% | %s |\n",
              pt.range.resolution_3db, pt.expected_range_3db, pt.range_error_pct, mark(range_ok));
  std::printf("| Azimuth -3 dB width | %.4f m | %.4f m | %+.2f%% | within 5%% | %s |\n",
              pt.azimuth.resolution_3db, pt.expected_azimuth_3db, pt.azimuth_error_pct,
              mark(az_ok));
  std::printf("| Range PSLR | %.2f dB | %.2f dB | %+.2f dB | within 1 dB | %s |\n",
              pt.range.pslr_db, pt.expected_pslr_db, pt.range.pslr_db - pt.expected_pslr_db,
              mark(pslr_ok));
  std::printf("| Azimuth PSLR | %.2f dB | %.2f dB | %+.2f dB | within 1 dB | %s |\n",
              pt.azimuth.pslr_db, pt.expected_pslr_db, pt.azimuth.pslr_db - pt.expected_pslr_db,
              mark(pslr_ok));
  std::printf("| Range ISLR | %.2f dB | -9.7 dB (infinite sinc) | | reported | |\n",
              pt.range.islr_db);
  std::printf("| Azimuth ISLR | %.2f dB | -9.7 dB (infinite sinc) | | reported | |\n\n",
              pt.azimuth.islr_db);

  const SpeckleReport sp = run_speckle_check(radar, geom, 80.0, 10.0, 1, threads);
  const bool enl_ok = std::abs(sp.stats.enl - 1.0) < 0.15;
  const bool ks_ok = sp.stats.exponential_at_1pct;
  all = all && enl_ok && ks_ok;
  std::printf("### Speckle (%zu scatterers, %zu pulses, %zu samples)\n\n", sp.scatterers, sp.pulses,
              sp.stats.n);
  std::printf("| Metric | Measured | Criterion | Result |\n|---|---|---|---|\n");
  std::printf("| ENL (mean^2 / var of intensity) | %.3f | 1 +/- 0.15 | %s |\n", sp.stats.enl,
              mark(enl_ok));
  std::printf(
      "| KS distance to exponential | D = %.4f, D* = %.3f | D* < 1.308 (1%% level) | %s "
      "|\n\n",
      sp.stats.ks_d, sp.stats.ks_modified, mark(ks_ok));

  std::printf("### Layover and shadow vs. closed form (incidence 35 deg)\n\n");
  std::printf(
      "| West slope | East slope | Layover expected | Layover bins | Shadow expected | "
      "Shadow bins | Result |\n|---|---|---|---|---|---|---|\n");
  for (const LayoverShadowCase& c : run_layover_shadow_checks(threads)) {
    all = all && c.pass();
    std::printf("| %.0f deg | %.0f deg | %s | %.2f%% | %s | %.2f%% | %s |\n", c.west_slope_deg,
                c.east_slope_deg, c.expect_layover ? "yes" : "no", 100.0 * c.layover_fraction,
                c.expect_shadow ? "yes" : "no", 100.0 * c.shadow_fraction, mark(c.pass()));
  }
  std::printf("\n");

  const DihedralReport dh = run_dihedral_check(threads);
  const bool dh_ok = std::abs(dh.measured_range - dh.predicted_range) <= dh.range_spacing &&
                     dh.double_to_single_db > 3.0;
  all = all && dh_ok;
  std::printf("### Dihedral (40 m box, 30 m tall)\n\n");
  std::printf(
      "| Corner range (theory) | Double-bounce peak | Error | Bin | Double vs. single "
      "bounce | Result |\n|---|---|---|---|---|---|\n");
  std::printf("| %.2f m | %.2f m | %+.2f m | %.1f m | %+.1f dB | %s |\n\n", dh.predicted_range,
              dh.measured_range, dh.measured_range - dh.predicted_range, dh.range_spacing,
              dh.double_to_single_db, mark(dh_ok));

  const DeterminismReport det = run_determinism_check(threads);
  all = all && det.geometric_identical && det.coherent_identical;
  std::printf("### Determinism (1 vs. %u threads)\n\n", det.threads_compared);
  std::printf("| Mode | Bit-identical | Result |\n|---|---|---|\n");
  std::printf("| Geometric | %s | %s |\n", det.geometric_identical ? "yes" : "no",
              mark(det.geometric_identical));
  std::printf("| Coherent | %s | %s |\n\n", det.coherent_identical ? "yes" : "no",
              mark(det.coherent_identical));

  std::printf("Overall: %s\n", all ? "all checks pass" : "SOME CHECKS FAILED");
  return all ? 0 : 1;
}

int run_irf(const cli::Args& args) {
  if (args.positional().size() != 2) {
    throw std::invalid_argument("irf needs one .npy file");
  }
  const double dx = args.num("dx", 1.0);
  const double dy = args.num("dy", 1.0);
  args.reject_unknown();
  const NpyArray a = read_npy(args.positional()[1]);
  if (a.shape.size() != 2) {
    throw std::invalid_argument("irf expects a 2-D array");
  }
  const std::size_t h = a.shape[0];
  const std::size_t w = a.shape[1];
  const auto peak = static_cast<std::size_t>(std::max_element(a.values.begin(), a.values.end()) -
                                             a.values.begin());
  std::vector<float> row(w);
  std::vector<float> col(h);
  for (std::size_t i = 0; i < w; ++i) {
    row[i] = static_cast<float>(a.values[(peak / w) * w + i]);
  }
  for (std::size_t j = 0; j < h; ++j) {
    col[j] = static_cast<float>(a.values[j * w + peak % w]);
  }
  for (const auto& [name, cut, spacing] :
       {std::tuple{"x (columns)", &row, dx}, std::tuple{"y (rows)", &col, dy}}) {
    const IrfMetrics m = analyze_irf(*cut, spacing);
    std::printf("%s: -3 dB width %.4f, PSLR %.2f dB, ISLR %.2f dB%s\n", name, m.resolution_3db,
                m.pslr_db, m.islr_db, m.valid ? "" : " (cut too short)");
  }
  return 0;
}

int run_speckle(const cli::Args& args) {
  if (args.positional().size() != 2) {
    throw std::invalid_argument("speckle needs one .npy file");
  }
  const bool is_intensity = args.has("intensity");
  args.reject_unknown();
  NpyArray a = read_npy(args.positional()[1]);
  if (!is_intensity) {
    for (double& v : a.values) {
      v *= v;
    }
  }
  const SpeckleStats s = speckle_statistics(a.values);
  std::printf("n %zu  mean %.6g  ENL %.4f  KS D %.4f  D* %.3f  exponential@5%% %s  @1%% %s\n", s.n,
              s.mean, s.enl, s.ks_d, s.ks_modified, s.exponential_at_5pct ? "yes" : "no",
              s.exponential_at_1pct ? "yes" : "no");
  return 0;
}

int run_iou(const cli::Args& args) {
  if (args.positional().size() != 3) {
    throw std::invalid_argument("iou needs two .npy files");
  }
  args.reject_unknown();
  const NpyArray a = read_npy(args.positional()[1]);
  const NpyArray b = read_npy(args.positional()[2]);
  std::vector<std::uint8_t> ma(a.values.size());
  std::vector<std::uint8_t> mb(b.values.size());
  std::transform(a.values.begin(), a.values.end(), ma.begin(),
                 [](double v) { return v != 0.0 ? 1 : 0; });
  std::transform(b.values.begin(), b.values.end(), mb.begin(),
                 [](double v) { return v != 0.0 ? 1 : 0; });
  std::printf("IoU %.4f\n", mask_iou(ma, mb));
  return 0;
}

}  // namespace

int main(int argc, char** argv) try {
  const cli::Args args(argc, argv, {"intensity", "help"});
  if (args.has("help") || args.positional().empty()) {
    std::cout << kUsage;
    return args.has("help") ? 0 : 1;
  }
  const std::string cmd = args.positional()[0];
  if (cmd == "theory") {
    return run_theory(args);
  }
  if (cmd == "irf") {
    return run_irf(args);
  }
  if (cmd == "speckle") {
    return run_speckle(args);
  }
  if (cmd == "iou") {
    return run_iou(args);
  }
  throw std::invalid_argument("unknown command '" + cmd + "'");
} catch (const std::exception& e) {
  std::fprintf(stderr, "bsar_eval: %s\n\n%s", e.what(), kUsage);
  return 2;
}
