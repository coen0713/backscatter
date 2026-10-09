#include "backscatter/sensor/radar.hpp"

#include <cmath>
#include <stdexcept>

namespace bsar {

std::vector<cdouble> make_chirp(const RadarParams& radar) {
  if (radar.sample_rate <= radar.bandwidth) {
    throw std::invalid_argument("make_chirp: sample rate must exceed the chirp bandwidth");
  }
  // Round before ceil so 10 us at 40 MHz is 400 samples, not 401.
  const auto n = static_cast<std::size_t>(
      std::ceil(std::round(radar.pulse_duration * radar.sample_rate * 1e6) * 1e-6));
  const double k = radar.chirp_rate();
  const double half = 0.5 * radar.pulse_duration;
  std::vector<cdouble> chirp(n);
  for (std::size_t i = 0; i < n; ++i) {
    const double t = static_cast<double>(i) / radar.sample_rate - half;
    chirp[i] = std::polar(1.0, kPi * k * t * t);
  }
  return chirp;
}

}  // namespace bsar
