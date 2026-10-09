#include "backscatter/eval/metrics.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace bsar {

IrfMetrics analyze_irf(std::span<const float> amp, double spacing) {
  IrfMetrics m;
  const std::size_t n = amp.size();
  if (n < 5) {
    return m;
  }
  const auto peak_it = std::max_element(amp.begin(), amp.end());
  const auto pk = static_cast<std::size_t>(peak_it - amp.begin());
  const double peak = *peak_it;
  if (!(peak > 0.0) || pk == 0 || pk + 1 >= n) {
    return m;
  }
  auto power = [&](std::size_t i) { return static_cast<double>(amp[i]) * amp[i]; };

  // Sub-sample peak position from a parabola through log-free power samples.
  const double pl = power(pk - 1);
  const double pc = power(pk);
  const double pr = power(pk + 1);
  const double denom = pl - 2 * pc + pr;
  m.peak_position = static_cast<double>(pk) + (denom != 0.0 ? 0.5 * (pl - pr) / denom : 0.0);
  m.peak_amplitude = peak;

  // Half-power crossings, linearly interpolated in power.
  const double half = 0.5 * pc;
  std::size_t l = pk;
  while (l > 0 && power(l) >= half) {
    --l;
  }
  std::size_t r = pk;
  while (r + 1 < n && power(r) >= half) {
    ++r;
  }
  if (power(l) >= half || power(r) >= half) {
    return m;  // cut too short to contain the main lobe
  }
  const double xl = static_cast<double>(l) + (half - power(l)) / (power(l + 1) - power(l));
  const double xr = static_cast<double>(r) - (half - power(r)) / (power(r - 1) - power(r));
  m.resolution_3db = (xr - xl) * spacing;

  // Main lobe = out to the first local minimum on each side.
  std::size_t nl = pk;
  while (nl > 0 && amp[nl - 1] < amp[nl]) {
    --nl;
  }
  std::size_t nr = pk;
  while (nr + 1 < n && amp[nr + 1] < amp[nr]) {
    ++nr;
  }
  double main_energy = 0.0;
  double side_energy = 0.0;
  double side_peak = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    const double p = power(i);
    if (i > nl && i < nr) {
      main_energy += p;
    } else {
      side_energy += p;
      side_peak = std::max(side_peak, p);
    }
  }
  m.pslr_db = 10.0 * std::log10(std::max(side_peak, 1e-300) / pc);
  m.islr_db = 10.0 * std::log10(std::max(side_energy, 1e-300) / main_energy);
  m.valid = nl > 0 && nr + 1 < n;
  return m;
}

SpeckleStats speckle_statistics(std::span<const double> intensity) {
  SpeckleStats s;
  s.n = intensity.size();
  if (s.n < 2) {
    return s;
  }
  double sum = 0.0;
  for (const double v : intensity) {
    sum += v;
  }
  s.mean = sum / static_cast<double>(s.n);
  double ss = 0.0;
  for (const double v : intensity) {
    ss += (v - s.mean) * (v - s.mean);
  }
  s.variance = ss / static_cast<double>(s.n - 1);
  s.enl = s.variance > 0.0 ? s.mean * s.mean / s.variance : 0.0;

  std::vector<double> sorted(intensity.begin(), intensity.end());
  std::sort(sorted.begin(), sorted.end());
  const auto nd = static_cast<double>(s.n);
  double d = 0.0;
  for (std::size_t i = 0; i < s.n; ++i) {
    const double cdf = 1.0 - std::exp(-sorted[i] / s.mean);
    const double lo = static_cast<double>(i) / nd;
    const double hi = static_cast<double>(i + 1) / nd;
    d = std::max({d, cdf - lo, hi - cdf});
  }
  s.ks_d = d;
  // Stephens (1974), exponential distribution with estimated mean.
  s.ks_modified = (d - 0.2 / nd) * (std::sqrt(nd) + 0.26 + 0.5 / std::sqrt(nd));
  s.exponential_at_5pct = s.ks_modified < 1.094;
  s.exponential_at_1pct = s.ks_modified < 1.308;
  return s;
}

double mask_iou(std::span<const std::uint8_t> a, std::span<const std::uint8_t> b) {
  if (a.size() != b.size()) {
    throw std::invalid_argument("mask_iou: masks differ in size");
  }
  std::size_t inter = 0;
  std::size_t uni = 0;
  for (std::size_t i = 0; i < a.size(); ++i) {
    const bool x = a[i] != 0;
    const bool y = b[i] != 0;
    inter += (x && y) ? 1 : 0;
    uni += (x || y) ? 1 : 0;
  }
  return uni == 0 ? 1.0 : static_cast<double>(inter) / static_cast<double>(uni);
}

}  // namespace bsar
