#include "backscatter/image/backprojection.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "backscatter/math/constants.hpp"
#include "backscatter/util/cpu.hpp"
#include "backscatter/util/parallel.hpp"

#include "backprojection_simd.hpp"

namespace bsar {
namespace {

using cacc = std::complex<double>;

/// Contribution of one pulse to one pixel.
inline cacc pulse_contribution(const CompressedData& d, const cfloat* row, const Vec3d& pulse,
                               const Vec3d& pixel, double inv_spacing, double phase_scale) {
  const double r = length(pixel - pulse);
  const double s = (r - d.near_range) * inv_spacing;
  if (!(s >= 0.0) || s >= static_cast<double>(d.num_samples - 1)) {
    return {};
  }
  const auto i = static_cast<std::size_t>(s);
  const auto f = static_cast<float>(s - static_cast<double>(i));
  const cfloat v = row[i] * (1.0f - f) + row[i + 1] * f;
  return cacc(v) * std::polar(1.0, phase_scale * r);
}

}  // namespace

BackprojectionKernel resolve_kernel(BackprojectionKernel requested) {
  if (requested == BackprojectionKernel::Auto || requested == BackprojectionKernel::Simd) {
    return cpu_has_avx2_fma() ? BackprojectionKernel::Simd : BackprojectionKernel::Blocked;
  }
  return requested;
}

Image<std::complex<float>> backproject(const CompressedData& data, const PixelGrid& grid,
                                       const BackprojectionConfig& config) {
  if (!grid.heights.empty() && grid.heights.size() != grid.width * grid.height) {
    throw std::invalid_argument("backproject: heights must have width * height entries");
  }
  Image<std::complex<float>> out(grid.width, grid.height);
  if (data.num_samples < 2 || data.num_pulses == 0) {
    return out;
  }
  const double inv_spacing = 1.0 / data.range_spacing;
  const double phase_scale = 4.0 * kPi / data.wavelength;
  const std::size_t n_pix = grid.width * grid.height;

  const BackprojectionKernel kernel = resolve_kernel(config.kernel);
  if (kernel == BackprojectionKernel::Naive) {
    parallel_for(0, grid.height, config.threads, [&](std::size_t j) {
      for (std::size_t i = 0; i < grid.width; ++i) {
        const Vec3d pixel = grid.position(i, j);
        cacc acc{};
        for (std::size_t k = 0; k < data.num_pulses; ++k) {
          acc += pulse_contribution(data, data.samples.data() + k * data.num_samples,
                                    data.positions[k], pixel, inv_spacing, phase_scale);
        }
        out.at(i, j) = std::complex<float>(acc);
      }
    });
    return out;
  }

  const std::size_t pix_block = std::max<std::size_t>(config.pixel_block, 1);
  const std::size_t pulse_block = std::max<std::size_t>(config.pulse_block, 1);
  const std::size_t n_tiles = (n_pix + pix_block - 1) / pix_block;

  if (kernel == BackprojectionKernel::Simd) {
    std::vector<double> pulse_xyz(3 * data.num_pulses);
    for (std::size_t k = 0; k < data.num_pulses; ++k) {
      pulse_xyz[3 * k] = data.positions[k].x;
      pulse_xyz[3 * k + 1] = data.positions[k].y;
      pulse_xyz[3 * k + 2] = data.positions[k].z;
    }
    const detail::BackprojectionInputs in{data.samples.data(), data.num_samples, pulse_xyz.data(),
                                          data.near_range,     inv_spacing,      phase_scale};
    parallel_for(0, n_tiles, config.threads, [&](std::size_t tile) {
      const std::size_t p0 = tile * pix_block;
      const std::size_t n = std::min(n_pix, p0 + pix_block) - p0;
      std::vector<double> px(n);
      std::vector<double> py(n);
      std::vector<double> pz(n);
      for (std::size_t p = 0; p < n; ++p) {
        const Vec3d v = grid.position((p0 + p) % grid.width, (p0 + p) / grid.width);
        px[p] = v.x;
        py[p] = v.y;
        pz[p] = v.z;
      }
      std::vector<double> acc(2 * n, 0.0);
      for (std::size_t k0 = 0; k0 < data.num_pulses; k0 += pulse_block) {
        const std::size_t k1 = std::min(data.num_pulses, k0 + pulse_block);
        const std::size_t done =
            detail::accumulate_avx2(in, px.data(), py.data(), pz.data(), n, k0, k1, acc.data());
        for (std::size_t p = done; p < n; ++p) {  // remainder pixels
          cacc a(acc[2 * p], acc[2 * p + 1]);
          for (std::size_t k = k0; k < k1; ++k) {
            a += pulse_contribution(data, data.samples.data() + k * data.num_samples,
                                    data.positions[k], {px[p], py[p], pz[p]}, inv_spacing,
                                    phase_scale);
          }
          acc[2 * p] = a.real();
          acc[2 * p + 1] = a.imag();
        }
      }
      for (std::size_t p = 0; p < n; ++p) {
        out.data[p0 + p] = std::complex<float>(cacc(acc[2 * p], acc[2 * p + 1]));
      }
    });
    return out;
  }

  parallel_for(0, n_tiles, config.threads, [&](std::size_t tile) {
    const std::size_t p0 = tile * pix_block;
    const std::size_t p1 = std::min(n_pix, p0 + pix_block);
    std::vector<Vec3d> pixels(p1 - p0);
    for (std::size_t p = p0; p < p1; ++p) {
      pixels[p - p0] = grid.position(p % grid.width, p / grid.width);
    }
    std::vector<cacc> acc(p1 - p0);
    for (std::size_t k0 = 0; k0 < data.num_pulses; k0 += pulse_block) {
      const std::size_t k1 = std::min(data.num_pulses, k0 + pulse_block);
      for (std::size_t p = 0; p < pixels.size(); ++p) {
        cacc a = acc[p];
        for (std::size_t k = k0; k < k1; ++k) {
          a += pulse_contribution(data, data.samples.data() + k * data.num_samples,
                                  data.positions[k], pixels[p], inv_spacing, phase_scale);
        }
        acc[p] = a;
      }
    }
    for (std::size_t p = p0; p < p1; ++p) {
      out.data[p] = std::complex<float>(acc[p - p0]);
    }
  });
  return out;
}

}  // namespace bsar
