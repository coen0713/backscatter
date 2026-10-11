#pragma once

#include <cstddef>
#include <cstdint>
#include <span>

#include "backscatter/image/image.hpp"
#include "backscatter/math/vec3.hpp"
#include "backscatter/scene/scene.hpp"
#include "backscatter/sensor/trajectory.hpp"
#include "backscatter/trace/geometric.hpp"

namespace bsar {

/// Slant-range products resampled onto ground points.
struct GeocodedProducts {
  Image<float> beta0;           // bilinear; 0 outside the image
  Image<std::uint8_t> layover;  // slant-range layover mask at the point's bin (nearest)
  Image<std::uint8_t> shadow;   // ground shadow: the sensor cannot see the point
  Image<std::uint8_t> valid;    // 1 where the point maps inside the image
  Image<float> range;           // zero-Doppler slant range [m] (diagnostics)
};

/// Range-Doppler geocoding: for each ground point (scene ENU, row-major
/// width x height; NaN marks no-data), solve for its zero-Doppler time and
/// slant range on `trajectory` and sample `image` there. This is the same
/// mapping SNAP's Range-Doppler Terrain Correction applies to real data, so
/// simulated and real products can be compared on one grid.
///
/// Masks follow the ground-geometry convention of SNAP's layover-shadow mask:
/// a point is in layover if its slant-range bin receives returns from more
/// than one part of the surface, and in shadow if it is occluded from the
/// sensor at its zero-Doppler time (an occlusion ray through `scene`).
GeocodedProducts geocode(const GeometricImage& image, const Trajectory& trajectory,
                         const Scene& scene, std::span<const Vec3d> points, std::size_t width,
                         std::size_t height, unsigned threads = 0);

}  // namespace bsar
