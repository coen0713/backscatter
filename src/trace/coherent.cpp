#include "backscatter/trace/coherent.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

#include "backscatter/math/constants.hpp"
#include "backscatter/math/rng.hpp"
#include "backscatter/trace/geometric.hpp"
#include "backscatter/util/parallel.hpp"

namespace bsar {
namespace {

// Fractional-delay interpolation kernel: Kaiser-windowed sinc, tabulated at
// kPhases sub-sample offsets.
constexpr int kHalfWidth = 16;
constexpr int kTaps = 2 * kHalfWidth;
constexpr int kPhases = 1024;
constexpr double kKaiserBeta = 6.0;

double bessel_i0(double x) {
  double sum = 1.0;
  double term = 1.0;
  for (int k = 1; k < 50; ++k) {
    term *= (x / (2.0 * k)) * (x / (2.0 * k));
    sum += term;
    if (term < 1e-17 * sum) {
      break;
    }
  }
  return sum;
}

const std::vector<double>& delay_kernel() {
  static const std::vector<double> table = [] {
    std::vector<double> t(static_cast<std::size_t>(kPhases + 1) * kTaps);
    const double norm = 1.0 / bessel_i0(kKaiserBeta);
    for (int ph = 0; ph <= kPhases; ++ph) {
      const double frac = static_cast<double>(ph) / kPhases;
      for (int j = 0; j < kTaps; ++j) {
        const double x = static_cast<double>(j - (kHalfWidth - 1)) - frac;  // n - position
        const double sinc = x == 0.0 ? 1.0 : std::sin(kPi * x) / (kPi * x);
        const double r = x / kHalfWidth;
        const double win =
            r * r < 1.0 ? bessel_i0(kKaiserBeta * std::sqrt(1.0 - r * r)) * norm : 0.0;
        t[static_cast<std::size_t>(ph) * kTaps + static_cast<std::size_t>(j)] = sinc * win;
      }
    }
    return t;
  }();
  return table;
}

/// Two-way azimuth antenna weight for a target at azimuth angle asin(sin_psi)
/// off boresight (boresight = zero-Doppler direction).
double beam_weight(const RadarParams& radar, double sin_psi) {
  const double x = radar.antenna_length * sin_psi / radar.wavelength();
  if (radar.beam == BeamPattern::Uniform) {
    return std::abs(x) <= 0.5 ? 1.0 : 0.0;
  }
  if (std::abs(x) >= 1.0) {
    return 0.0;
  }
  const double s = x == 0.0 ? 1.0 : std::sin(kPi * x) / (kPi * x);
  return s * s;
}

}  // namespace

std::vector<Scatterer> sample_scatterers(const Scene& scene, const Vec3d& illuminator,
                                         const ScattererConfig& config,
                                         const ScatteringModel& model) {
  if (!scene.built()) {
    throw std::invalid_argument("sample_scatterers: scene is not built");
  }
  const CounterRng rng(config.seed);
  const std::size_t n_tri = scene.mesh.num_triangles();
  constexpr std::size_t kChunk = 1024;
  const std::size_t n_chunks = (n_tri + kChunk - 1) / kChunk;
  std::vector<std::vector<Scatterer>> per_chunk(n_chunks);
  const double eps = 1e-3 + 1e-6 * length(Vec3d(scene.bounds().extent()));

  parallel_for(0, n_chunks, 0, [&](std::size_t chunk) {
    std::vector<Scatterer>& out = per_chunk[chunk];
    const std::size_t end = std::min(n_tri, (chunk + 1) * kChunk);
    for (std::size_t t = chunk * kChunk; t < end; ++t) {
      const auto [fa, fb, fc] = scene.mesh.triangle(t);
      const Vec3d a(fa);
      const Vec3d b(fb);
      const Vec3d c(fc);
      const Vec3d cr = cross(b - a, c - a);
      const double area = 0.5 * length(cr);
      if (area <= 0.0) {
        continue;
      }
      const Vec3d n = normalize(cr);
      const double expected = config.density * area;
      const double whole = std::floor(expected);
      const double u_round = rng.uniform4(t, ~std::uint64_t{0})[0];
      const auto count = static_cast<std::size_t>(whole) + (u_round < expected - whole ? 1u : 0u);
      if (count == 0) {
        continue;
      }
      const Material& mat = scene.material_of(static_cast<std::uint32_t>(t));
      for (std::size_t k = 0; k < count; ++k) {
        const auto u = rng.uniform4(t, k);
        const double r1 = std::sqrt(u[0]);
        const Vec3d pos = a * (1.0 - r1) + b * (r1 * (1.0 - u[1])) + c * (r1 * u[1]);
        const Vec3d to_ill = illuminator - pos;
        const Vec3d s = normalize(to_ill);
        const double cos_inc = dot(n, s);
        if (cos_inc <= 0.0) {
          continue;  // facets are one-sided for coherent scattering
        }
        if (config.check_visibility) {
          const Ray shadow_ray{Vec3f(pos + n * eps), Vec3f(s), 0.0f,
                               std::numeric_limits<float>::infinity()};
          if (scene.occluded(shadow_ray)) {
            continue;
          }
        }
        const double sample_area = area / static_cast<double>(count);
        const double sigma0 = model.backscatter(mat, cos_inc, config.polarization);
        out.emplace_back(pos, std::sqrt(sigma0 * sample_area));

        // Follow the specular reflection of the illuminating ray. Weighting
        // matches the geometric integrator: the beam cross-section carried
        // by this sample is sample_area * cos_inc, spread over
        // 1 / cos(incidence) on each surface it lands on.
        const double tube = sample_area * cos_inc;
        Vec3d dir = -s;
        Vec3d cur = pos;
        Vec3d n_cur = n;
        const Material* m_cur = &mat;
        double cos_cur = cos_inc;
        double weight = 1.0;
        double internal = 0.0;
        for (int bounce = 2; bounce <= config.max_bounces; ++bounce) {
          weight *= model.specular_reflectance(*m_cur, cos_cur, config.polarization);
          if (weight < 1e-6) {
            break;
          }
          dir = reflect(dir, n_cur);
          const Ray ray{Vec3f(cur + n_cur * eps), Vec3f(dir), 0.0f,
                        std::numeric_limits<float>::infinity()};
          Hit hit;
          if (!scene.intersect(ray, hit)) {
            break;
          }
          const Vec3d next = Vec3d(ray.origin) + Vec3d(ray.dir) * static_cast<double>(hit.t);
          internal += length(next - cur);
          Vec3d nn(scene.mesh.normal(hit.prim));
          if (dot(nn, dir) > 0.0) {
            nn = -nn;
          }
          const Material& mn = scene.material_of(hit.prim);
          const double cos_in = std::max(-dot(nn, dir), 1e-6);
          const Vec3d back = normalize(illuminator - next);
          if (dot(nn, back) > 0.0) {
            const Ray shadow_ray{Vec3f(next + nn * eps), Vec3f(back), 0.0f,
                                 std::numeric_limits<float>::infinity()};
            if (!config.check_visibility || !scene.occluded(shadow_ray)) {
              const double rcs =
                  weight * model.bistatic(mn, nn, -dir, back, config.polarization) * tube / cos_in;
              if (rcs > 0.0) {
                out.emplace_back(pos, next, internal, std::sqrt(rcs), bounce);
              }
            }
          }
          cur = next;
          n_cur = nn;
          m_cur = &mn;
          cos_cur = cos_in;
        }
      }
    }
  });

  std::vector<Scatterer> all;
  for (auto& v : per_chunk) {
    all.insert(all.end(), v.begin(), v.end());
  }
  return all;
}

EchoWindow auto_echo_window(const Aabb& bounds, const Trajectory& trajectory,
                            const RadarParams& radar) {
  double t_lo = std::numeric_limits<double>::infinity();
  double t_hi = -t_lo;
  double r_lo = std::numeric_limits<double>::infinity();
  double r_hi = -r_lo;
  for (int c = 0; c < 8; ++c) {
    const Vec3d corner(bounds.corner(c));
    const double t = find_zero_doppler_time(trajectory, corner);
    const double r = length(trajectory.position(t) - corner);
    t_lo = std::min(t_lo, t);
    t_hi = std::max(t_hi, t);
    r_lo = std::min(r_lo, r);
    r_hi = std::max(r_hi, r);
  }
  // The main lobe of the two-way pattern spans +/- lambda / L for the sinc^2
  // pattern and +/- lambda / (2L) for the uniform one.
  const double half_angle =
      (radar.beam == BeamPattern::Uniform ? 0.5 : 1.0) * radar.wavelength() / radar.antenna_length;
  const double half_aperture = r_hi * std::tan(half_angle);
  const double speed = length(trajectory.velocity(0.5 * (t_lo + t_hi)));
  const double margin = 4.0 * kSpeedOfLight / (2.0 * radar.bandwidth);
  EchoWindow w;
  w.t_start = t_lo - half_aperture / speed;
  w.t_end = t_hi + half_aperture / speed;
  w.near_range = r_lo - margin;
  w.far_range = std::hypot(r_hi, half_aperture) + margin;
  return w;
}

RawData synthesize_raw(std::span<const Scatterer> scatterers, const Trajectory& trajectory,
                       const RadarParams& radar, const EchoWindow& window, unsigned threads) {
  if (!(window.far_range > window.near_range) || !(window.t_end >= window.t_start)) {
    throw std::invalid_argument("synthesize_raw: empty echo window");
  }
  const double fs = radar.sample_rate;
  const double lambda = radar.wavelength();
  const std::vector<cdouble> chirp = make_chirp(radar);
  const std::size_t nc = chirp.size();
  const auto n_imp = static_cast<std::size_t>(std::ceil((window.far_range - window.near_range) *
                                                        2.0 / kSpeedOfLight * fs)) +
                     1;
  const std::size_t n_raw = n_imp + nc - 1;
  const std::size_t m = next_pow2(n_imp + nc);
  const double tau_near = 2.0 * window.near_range / kSpeedOfLight;

  std::vector<cdouble> chirp_spec(m);
  std::copy(chirp.begin(), chirp.end(), chirp_spec.begin());
  fft_inplace(chirp_spec);

  RawData raw;
  raw.radar = radar;
  raw.near_range = window.near_range;
  raw.window_samples = n_imp;
  raw.num_samples = n_raw;
  raw.num_pulses =
      static_cast<std::size_t>(std::floor((window.t_end - window.t_start) * radar.prf)) + 1;
  raw.times.resize(raw.num_pulses);
  raw.positions.resize(raw.num_pulses);
  for (std::size_t k = 0; k < raw.num_pulses; ++k) {
    raw.times[k] = window.t_start + static_cast<double>(k) / radar.prf;
    raw.positions[k] = trajectory.position(raw.times[k]);
  }
  raw.samples.assign(raw.num_pulses * n_raw, cfloat{});

  const std::vector<double>& kernel = delay_kernel();
  const double phase_scale = -4.0 * kPi / lambda;

  parallel_for(0, raw.num_pulses, threads, [&](std::size_t k) {
    const Vec3d p = raw.positions[k];
    const Vec3d ev = normalize(trajectory.velocity(raw.times[k]));
    std::vector<cdouble> buf(m, cdouble{});
    for (const Scatterer& s : scatterers) {
      const Vec3d d = (s.bounces == 1 ? s.position : (s.position + s.exit) * 0.5) - p;
      const double r = s.bounces == 1 ? length(d) : s.range_from(p);
      const double w = beam_weight(radar, dot(d, ev) / length(d));
      if (w == 0.0) {
        continue;
      }
      const double x = (2.0 * r / kSpeedOfLight - tau_near) * fs;
      if (x < -kHalfWidth || x > static_cast<double>(n_imp + kHalfWidth)) {
        continue;
      }
      const cdouble amp = std::polar(s.amplitude * w, phase_scale * r);
      auto i0 = static_cast<std::ptrdiff_t>(std::floor(x));
      auto ph = static_cast<std::ptrdiff_t>(std::lround((x - static_cast<double>(i0)) * kPhases));
      const double* taps = kernel.data() + ph * kTaps;
      for (std::ptrdiff_t j = 0; j < kTaps; ++j) {
        const std::ptrdiff_t n = i0 + j - (kHalfWidth - 1);
        if (n >= 0 && n < static_cast<std::ptrdiff_t>(n_imp)) {
          buf[static_cast<std::size_t>(n)] += amp * taps[j];
        }
      }
    }
    fft_inplace(buf);
    for (std::size_t i = 0; i < m; ++i) {
      buf[i] *= chirp_spec[i];
    }
    fft_inplace(buf, true);
    cfloat* out = raw.samples.data() + k * n_raw;
    for (std::size_t i = 0; i < n_raw; ++i) {
      out[i] = cfloat(buf[i]);
    }
  });
  return raw;
}

CompressedData range_compress(const RawData& raw, int upsample, unsigned threads) {
  if (upsample < 1) {
    throw std::invalid_argument("range_compress: upsample must be >= 1");
  }
  const std::vector<cdouble> chirp = make_chirp(raw.radar);
  const std::size_t nc = chirp.size();
  const std::size_t m = next_pow2(raw.num_samples + nc);
  const auto up = static_cast<std::size_t>(upsample);
  const std::size_t mu = m * up;

  std::vector<cdouble> ref(m);
  std::copy(chirp.begin(), chirp.end(), ref.begin());
  fft_inplace(ref);
  for (auto& v : ref) {
    v = std::conj(v);
  }

  CompressedData out;
  out.num_pulses = raw.num_pulses;
  out.num_samples = raw.window_samples * up;
  out.positions = raw.positions;
  out.near_range = raw.near_range;
  out.range_spacing = kSpeedOfLight / (2.0 * raw.radar.sample_rate * static_cast<double>(up));
  out.wavelength = raw.radar.wavelength();
  out.samples.assign(out.num_pulses * out.num_samples, cfloat{});
  const double scale = static_cast<double>(up) / static_cast<double>(nc);

  parallel_for(0, raw.num_pulses, threads, [&](std::size_t k) {
    std::vector<cdouble> spec(m, cdouble{});
    const cfloat* in = raw.samples.data() + k * raw.num_samples;
    for (std::size_t i = 0; i < raw.num_samples; ++i) {
      spec[i] = cdouble(in[i]);
    }
    fft_inplace(spec);
    // Matched filter, then zero-insert in the middle of the spectrum.
    std::vector<cdouble> wide(mu, cdouble{});
    const std::size_t half = m / 2;
    for (std::size_t i = 0; i < half; ++i) {
      wide[i] = spec[i] * ref[i];
      wide[mu - half + i] = spec[half + i] * ref[half + i];
    }
    fft_inplace(wide, true);
    cfloat* dst = out.samples.data() + k * out.num_samples;
    for (std::size_t i = 0; i < out.num_samples; ++i) {
      dst[i] = cfloat(wide[i] * scale);
    }
  });
  return out;
}

}  // namespace bsar
