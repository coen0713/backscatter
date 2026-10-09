#pragma once

#include "backscatter/image/image.hpp"
#include "backscatter/scene/scene.hpp"

namespace bsar {

struct HillshadeConfig {
  double sun_azimuth_deg = 315.0;   // clockwise from north, GDAL default
  double sun_elevation_deg = 45.0;  // GDAL default
  double pixel_size = 0.0;          // metres; 0 picks ~512 px across
  bool cast_shadows = false;        // GDAL's hillshade has no cast shadows
  unsigned threads = 0;
};

/// Optical sanity check for the ray tracing core: an orthographic top-down
/// render of cos(angle between surface normal and sun), i.e. the quantity
/// `gdaldem hillshade` computes (scaled to [0, 1]). Row 0 is the northern edge.
Image<float> render_hillshade(const Scene& scene, const HillshadeConfig& config);

}  // namespace bsar
