#pragma once

#include <complex>
#include <cstddef>
#include <vector>

#include "backscatter/image/image.hpp"
#include "backscatter/math/vec3.hpp"
#include "backscatter/trace/coherent.hpp"

namespace bsar {

/// Output pixel positions: pixel (i, j) sits at origin + i * step_x + j * step_y.
/// If `heights` is non-empty (width * height values) it replaces the z
/// coordinate, e.g. to focus directly onto a DEM.
struct PixelGrid {
  std::size_t width = 0;
  std::size_t height = 0;
  Vec3d origin;
  Vec3d step_x{1, 0, 0};
  Vec3d step_y{0, 1, 0};
  std::vector<double> heights;

  [[nodiscard]] Vec3d position(std::size_t i, std::size_t j) const {
    Vec3d p = origin + step_x * static_cast<double>(i) + step_y * static_cast<double>(j);
    if (!heights.empty()) {
      p.z = heights[j * width + i];
    }
    return p;
  }
};

enum class BackprojectionKernel {
  Naive,    // pixel-major: every pulse for one pixel, then the next pixel
  Blocked,  // tiles of pixels x blocks of pulses, to keep pulse data in cache
};

struct BackprojectionConfig {
  BackprojectionKernel kernel = BackprojectionKernel::Blocked;
  unsigned threads = 0;
  std::size_t pixel_block = 256;
  std::size_t pulse_block = 32;
};

/// Time-domain backprojection: for each pixel, sum the range-compressed
/// echoes of every pulse at the pixel's range, phase-corrected by
/// exp(+i 4 pi R / lambda). O(pixels x pulses). Every pixel accumulates its
/// pulses in index order in both kernels, so results are bit-identical
/// across kernels and thread counts.
Image<std::complex<float>> backproject(const CompressedData& data, const PixelGrid& grid,
                                       const BackprojectionConfig& config = {});

}  // namespace bsar
