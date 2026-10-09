#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include "backscatter/eval/metrics.hpp"
#include "backscatter/sensor/radar.hpp"

namespace bsar {

/// The self-contained validation suite of docs/validation.md section 1.
/// Each check is a function so the unit tests and `bsar_eval theory` share
/// one implementation (and the README numbers come from the same code).

/// Compact airborne C-band system used by the coherent checks: 100 MHz
/// chirp (1.5 m slant resolution), 2 m antenna (1 m azimuth resolution).
RadarParams theory_radar();

struct CoherentGeometry {
  double altitude = 3000.0;      // m
  double ground_range = 4000.0;  // m, from nadir track to scene centre
  double speed = 100.0;          // m/s
};

struct PointTargetReport {
  IrfMetrics range;
  IrfMetrics azimuth;
  double expected_range_3db = 0.0;    // 0.886 c / (2B)
  double expected_azimuth_3db = 0.0;  // 0.886 L / 2
  double range_error_pct = 0.0;
  double azimuth_error_pct = 0.0;
  double expected_pslr_db = -13.26;  // unweighted sinc
  std::size_t pulses = 0;
};

PointTargetReport run_point_target_check(const RadarParams& radar, const CoherentGeometry& geom,
                                         unsigned threads = 0);

struct SpeckleReport {
  SpeckleStats stats;
  std::size_t scatterers = 0;
  std::size_t pulses = 0;
};

/// Homogeneous flat patch, single look; intensity sampled every two
/// resolution cells in each direction so samples are close to independent.
SpeckleReport run_speckle_check(const RadarParams& radar, const CoherentGeometry& geom,
                                double patch_size = 80.0, double density = 10.0,
                                std::uint64_t seed = 1, unsigned threads = 0);

struct LayoverShadowCase {
  double incidence_deg = 0.0;
  double west_slope_deg = 0.0;  // faces the (west-side, right-looking) sensor
  double east_slope_deg = 0.0;  // faces away
  bool expect_layover = false;  // closed form: west slope > incidence
  bool expect_shadow = false;   // closed form: east slope > 90 deg - incidence
  double layover_fraction = 0.0;
  double shadow_fraction = 0.0;
  [[nodiscard]] bool pass() const {
    return (layover_fraction > 0.0) == expect_layover && (shadow_fraction > 0.0) == expect_shadow;
  }
};

std::vector<LayoverShadowCase> run_layover_shadow_checks(unsigned threads = 0);

struct DihedralReport {
  double predicted_range = 0.0;  // slant range of the wall-ground corner
  double measured_range = 0.0;   // peak of the double-bounce image
  double range_spacing = 0.0;
  double double_to_single_db = 0.0;  // double-bounce peak vs. brightest single-bounce bin
};

DihedralReport run_dihedral_check(unsigned threads = 0);

struct DeterminismReport {
  bool geometric_identical = false;
  bool coherent_identical = false;
  unsigned threads_compared = 0;
};

/// Same inputs on 1 thread vs. `threads`: outputs must match bit for bit.
DeterminismReport run_determinism_check(unsigned threads = 0);

}  // namespace bsar
