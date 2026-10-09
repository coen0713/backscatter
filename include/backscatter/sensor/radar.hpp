#pragma once

#include <vector>

#include "backscatter/math/constants.hpp"
#include "backscatter/math/fft.hpp"

namespace bsar {

enum class BeamPattern {
  Uniform,  // hard cut-off at +/- lambda / (2 L): uniform azimuth weighting
  Sinc2,    // two-way sinc^2 amplitude pattern of a uniform aperture, main lobe only
};

/// Radar system parameters for coherent simulation. Defaults are close to
/// Sentinel-1 stripmap (C band).
struct RadarParams {
  double carrier_frequency = 5.405e9;  // Hz
  double bandwidth = 50e6;             // Hz, chirp bandwidth
  double pulse_duration = 20e-6;       // s
  double sample_rate = 60e6;           // Hz, complex baseband sampling
  double prf = 1700.0;                 // Hz
  double antenna_length = 12.3;        // m, azimuth
  BeamPattern beam = BeamPattern::Uniform;

  [[nodiscard]] double wavelength() const { return kSpeedOfLight / carrier_frequency; }
  [[nodiscard]] double chirp_rate() const { return bandwidth / pulse_duration; }
};

/// Baseband linear FM chirp s(t) = exp(i pi K (t - T/2)^2), t in [0, T),
/// sampled at `sample_rate`. Sweeps -B/2 .. +B/2 with unit amplitude.
std::vector<cdouble> make_chirp(const RadarParams& radar);

}  // namespace bsar
