// Python bindings: backscatter as a synthetic SAR data generator returning
// NumPy arrays.

#include <complex>
#include <cstring>
#include <memory>
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include "backscatter/eval/theory_checks.hpp"
#include "backscatter/image/backprojection.hpp"
#include "backscatter/scene/heightmap.hpp"
#include "backscatter/scene/scene.hpp"
#include "backscatter/sensor/sar_geometry.hpp"
#include "backscatter/trace/coherent.hpp"
#include "backscatter/trace/geometric.hpp"
#include "backscatter/trace/hillshade.hpp"

namespace nb = nanobind;
using namespace bsar;

namespace {

/// Move an Image into a NumPy array that owns the data.
template <typename T>
nb::ndarray<nb::numpy, T, nb::ndim<2>> to_numpy(Image<T>&& img) {
  auto* owned = new std::vector<T>(std::move(img.data));
  nb::capsule deleter(owned, [](void* p) noexcept { delete static_cast<std::vector<T>*>(p); });
  return nb::ndarray<nb::numpy, T, nb::ndim<2>>(owned->data(), {img.height, img.width}, deleter);
}

Scene scene_from_heights(nb::ndarray<const float, nb::ndim<2>, nb::c_contig> heights,
                         double cell_size, const std::string& material) {
  Heightmap hm;
  hm.ny = heights.shape(0);
  hm.nx = heights.shape(1);
  hm.dx = hm.dy = cell_size;
  hm.x0 = -0.5 * cell_size * static_cast<double>(hm.nx - 1);
  hm.y0 = -0.5 * cell_size * static_cast<double>(hm.ny - 1);
  hm.z.resize(hm.nx * hm.ny);
  // Row 0 of the array is the northern edge (image convention).
  for (std::size_t r = 0; r < hm.ny; ++r) {
    std::memcpy(hm.z.data() + (hm.ny - 1 - r) * hm.nx, heights.data() + r * hm.nx,
                hm.nx * sizeof(float));
  }
  Scene scene;
  scene.mesh = heightmap_to_mesh(hm, scene.material_id(material));
  scene.build();
  return scene;
}

LinearTrajectory track_for(const Scene& scene, const std::string& platform, double incidence,
                           const std::string& look) {
  const bool airborne = platform == "airborne";
  return make_side_looking_track(scene.bounds(), airborne ? 6000.0 : 693e3,
                                 airborne ? 150.0 : 7600.0, incidence,
                                 look == "left" ? LookSide::Left : LookSide::Right);
}

nb::dict geometric_to_dict(GeometricImage&& img) {
  nb::dict d;
  d["near_range"] = img.near_range;
  d["range_spacing"] = img.range_spacing;
  d["t_start"] = img.t_start;
  d["line_interval"] = img.line_interval;
  d["intensity"] = to_numpy(std::move(img.intensity));
  d["layover"] = to_numpy(std::move(img.layover));
  d["shadow"] = to_numpy(std::move(img.shadow));
  nb::list bounces;
  for (auto& b : img.bounce) {
    bounces.append(to_numpy(std::move(b)));
  }
  d["bounce"] = bounces;
  return d;
}

GeometricConfig geometric_config(double range_spacing, double azimuth_spacing, int bounces,
                                 unsigned threads) {
  GeometricConfig cfg;
  cfg.range_spacing = range_spacing;
  cfg.azimuth_spacing = azimuth_spacing;
  cfg.max_bounces = bounces;
  cfg.threads = threads;
  return cfg;
}

}  // namespace

NB_MODULE(backscatter, m) {
  m.doc() = "backscatter: a physically based SAR ray tracer";

  m.def(
      "render_synthetic",
      [](const std::string& scene_name, double incidence, double range_spacing,
         double azimuth_spacing, int bounces, const std::string& platform, const std::string& look,
         unsigned threads) {
        const Scene scene = make_synthetic_scene(scene_name);
        const auto track = track_for(scene, platform, incidence, look);
        GeometricImage img;
        {
          const nb::gil_scoped_release release;
          img = render_geometric(
              scene, track, geometric_config(range_spacing, azimuth_spacing, bounces, threads));
        }
        return geometric_to_dict(std::move(img));
      },
      nb::arg("scene") = "city", nb::arg("incidence") = 35.0, nb::arg("range_spacing") = 5.0,
      nb::arg("azimuth_spacing") = 5.0, nb::arg("bounces") = 3, nb::arg("platform") = "spaceborne",
      nb::arg("look") = "right", nb::arg("threads") = 0,
      "Geometric SAR render of a named synthetic scene (plane, ridge, dihedral, city, "
      "mountains). Returns a dict of NumPy arrays (rows = azimuth, columns = slant range).");

  m.def(
      "render_heightmap",
      [](nb::ndarray<const float, nb::ndim<2>, nb::c_contig> heights, double cell_size,
         double incidence, double range_spacing, double azimuth_spacing, int bounces,
         const std::string& material, const std::string& platform, const std::string& look,
         unsigned threads) {
        const Scene scene = scene_from_heights(heights, cell_size, material);
        const auto track = track_for(scene, platform, incidence, look);
        GeometricImage img;
        {
          const nb::gil_scoped_release release;
          img = render_geometric(
              scene, track, geometric_config(range_spacing, azimuth_spacing, bounces, threads));
        }
        return geometric_to_dict(std::move(img));
      },
      nb::arg("heights"), nb::arg("cell_size"), nb::arg("incidence") = 35.0,
      nb::arg("range_spacing") = 5.0, nb::arg("azimuth_spacing") = 5.0, nb::arg("bounces") = 3,
      nb::arg("material") = "dry_soil", nb::arg("platform") = "spaceborne",
      nb::arg("look") = "right", nb::arg("threads") = 0,
      "Geometric SAR render of a float32 height grid (row 0 = north, metres).");

  m.def(
      "hillshade",
      [](nb::ndarray<const float, nb::ndim<2>, nb::c_contig> heights, double cell_size,
         double azimuth, double elevation) {
        const Scene scene = scene_from_heights(heights, cell_size, "dry_soil");
        HillshadeConfig cfg;
        cfg.sun_azimuth_deg = azimuth;
        cfg.sun_elevation_deg = elevation;
        cfg.pixel_size = cell_size;
        return to_numpy(render_hillshade(scene, cfg));
      },
      nb::arg("heights"), nb::arg("cell_size"), nb::arg("azimuth") = 315.0,
      nb::arg("elevation") = 45.0, "Ray-traced hillshade in [0, 1] (comparable to gdaldem).");

  m.def(
      "simulate_slc",
      [](nb::ndarray<const float, nb::ndim<2>, nb::c_contig> heights, double cell_size,
         double pixel_size, double density, std::uint64_t seed, double incidence,
         unsigned threads) {
        const Scene scene = scene_from_heights(heights, cell_size, "dry_soil");
        const auto track = track_for(scene, "airborne", incidence, "right");
        RadarParams radar = theory_radar();
        radar.prf = 4.0 * 150.0 / radar.antenna_length;
        ScattererConfig sc;
        sc.density = density;
        sc.seed = seed;
        Image<std::complex<float>> slc;
        {
          const nb::gil_scoped_release release;
          const double t_c = find_zero_doppler_time(track, Vec3d(scene.bounds().center()));
          const auto scatterers =
              sample_scatterers(scene, track.position(t_c), sc, FacetScatteringModel{});
          const EchoWindow w = auto_echo_window(scene.bounds(), track, radar);
          const CompressedData comp =
              range_compress(synthesize_raw(scatterers, track, radar, w, threads), 8, threads);
          const Aabb& b = scene.bounds();
          PixelGrid grid;
          grid.width = static_cast<std::size_t>(b.extent().x / pixel_size);
          grid.height = static_cast<std::size_t>(b.extent().y / pixel_size);
          grid.origin = {b.lo.x + 0.5 * pixel_size, b.hi.y - 0.5 * pixel_size, 0.0};
          grid.step_x = {pixel_size, 0.0, 0.0};
          grid.step_y = {0.0, -pixel_size, 0.0};
          grid.heights.resize(grid.width * grid.height);
          for (std::size_t j = 0; j < grid.height; ++j) {
            for (std::size_t i = 0; i < grid.width; ++i) {
              const Vec3d p = grid.position(i, j);
              Ray down{Vec3f(Vec3d{p.x, p.y, b.hi.z + 1.0}), {0.0f, 0.0f, -1.0f}};
              Hit hit;
              grid.heights[j * grid.width + i] =
                  scene.intersect(down, hit) ? b.hi.z + 1.0 - hit.t : static_cast<double>(b.lo.z);
            }
          }
          BackprojectionConfig bp;
          bp.threads = threads;
          slc = backproject(comp, grid, bp);
        }
        return to_numpy(std::move(slc));
      },
      nb::arg("heights"), nb::arg("cell_size"), nb::arg("pixel_size") = 0.5,
      nb::arg("density") = 4.0, nb::arg("seed") = 1, nb::arg("incidence") = 53.0,
      nb::arg("threads") = 0,
      "Coherent (airborne C-band) simulation of a height grid, focused by backprojection onto "
      "the surface. Returns a complex64 single-look complex image (row 0 = north).");

  m.def(
      "point_target_check",
      [](unsigned threads) {
        const PointTargetReport r =
            run_point_target_check(theory_radar(), CoherentGeometry{}, threads);
        nb::dict d;
        d["range_3db"] = r.range.resolution_3db;
        d["azimuth_3db"] = r.azimuth.resolution_3db;
        d["expected_range_3db"] = r.expected_range_3db;
        d["expected_azimuth_3db"] = r.expected_azimuth_3db;
        d["range_pslr_db"] = r.range.pslr_db;
        d["azimuth_pslr_db"] = r.azimuth.pslr_db;
        d["range_islr_db"] = r.range.islr_db;
        d["azimuth_islr_db"] = r.azimuth.islr_db;
        return d;
      },
      nb::arg("threads") = 0, "Run the point-target resolution/sidelobe check.");
}
