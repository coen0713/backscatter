// bsar_render: scene + acquisition geometry -> SAR image.

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "backscatter/image/backprojection.hpp"
#include "backscatter/image/geocode.hpp"
#include "backscatter/image/io.hpp"
#include "backscatter/scene/buildings.hpp"
#include "backscatter/scene/heightmap.hpp"
#include "backscatter/scene/scene.hpp"
#include "backscatter/sensor/orbit.hpp"
#include "backscatter/sensor/passes.hpp"
#include "backscatter/sensor/sar_geometry.hpp"
#include "backscatter/sensor/trajectory.hpp"
#include "backscatter/trace/coherent.hpp"
#include "backscatter/trace/geometric.hpp"
#include "backscatter/trace/hillshade.hpp"

#include "cli.hpp"

using namespace bsar;

namespace {

constexpr const char* kUsage = R"(usage: bsar_render [options]

Scene
  --scene NAME|dem:PATH   synthetic scene (plane, ridge, dihedral, city, mountains)
                          or a DEM: ESRI ASCII grid (.asc), or any GDAL raster
                          when built with BSAR_WITH_GDAL       [default: city]
  --geographic            the .asc DEM is in degrees (lon/lat), e.g. Copernicus
  --buildings PATH        building footprints (see data/fetch_osm_buildings.py)
  --classes PATH          per-sample class grid on the DEM's grid (.asc), with
  --material-table PATH   class -> material table (validation/materials_worldcover.txt)
  --origin LAT,LON,H      geodetic origin of the scene ENU frame (required for
                          geographic DEMs and orbit files)

Acquisition
  --orbit PATH            Sentinel-1 EOF orbit file (default: straight track)
  --time UTC              restrict the orbit to +/-10 min of this time
  --platform P            synthetic track: spaceborne (693 km, 7.6 km/s) or
                          airborne (6 km, 150 m/s)          [default: spaceborne]
  --incidence DEG         synthetic track incidence at scene centre [35]
  --look right|left       synthetic track look side [right]

Rendering
  --mode M                geometric | coherent | hillshade [geometric]
  --range-spacing M       slant-range bin [m]           [geometric: 5]
  --azimuth-spacing M     line spacing [m]              [geometric: 5]
  --rays-per-bin X        ray oversampling              [4]
  --bounces N             maximum bounces               [3]
  --pol HH|VV             polarisation                  [VV]
  --geocode               geometric + DEM: also write beta0, layover, shadow and
                          valid masks on the DEM's grid (_geo_*.asc, with the
                          DEM's .prj copied alongside)
  --density D             coherent: scatterers per m^2  [airborne 4, spaceborne 0.2]
  --pixel M               coherent / hillshade: output pixel size [m]
  --sun-azimuth DEG       hillshade: sun azimuth, clockwise from north [315]
  --sun-elevation DEG     hillshade: sun elevation [45]
  --seed N                coherent: scatterer seed      [1]
  --threads N             0 = all cores                 [0]
  --out PREFIX            output path prefix            [bsar]
)";

struct Outputs {
  std::string prefix;
  [[nodiscard]] std::string path(const std::string& suffix) const { return prefix + suffix; }
};

void write_overlay(const Outputs& out, const GeometricImage& img) {
  const auto gray = to_db_u8_auto(img.intensity, 35.0);
  std::vector<std::uint8_t> rgb(gray.size() * 3);
  // Radar brightness in gray; layover tinted red, shadow tinted blue.
  auto blend = [](int v, int tint) { return static_cast<std::uint8_t>((v * 3 + tint * 2) / 5); };
  for (std::size_t i = 0; i < gray.size(); ++i) {
    const int v = gray.data[i];
    int r = v;
    int g = v;
    int b = v;
    if (img.layover.data[i] != 0) {
      r = blend(v, 255);
      g = blend(v, 40);
      b = blend(v, 40);
    } else if (img.shadow.data[i] != 0) {
      r = blend(v, 20);
      g = blend(v, 60);
      b = blend(v, 200);
    }
    rgb[3 * i] = static_cast<std::uint8_t>(r);
    rgb[3 * i + 1] = static_cast<std::uint8_t>(g);
    rgb[3 * i + 2] = static_cast<std::uint8_t>(b);
  }
  write_png(out.path("_overlay.png"), gray.width, gray.height, 3, rgb);
}

double seconds_since(std::chrono::steady_clock::time_point t0) {
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

}  // namespace

int main(int argc, char** argv) try {
  const cli::Args args(argc, argv, {"geographic", "geocode", "help"});
  if (args.has("help")) {
    std::cout << kUsage;
    return 0;
  }
  const auto t0 = std::chrono::steady_clock::now();
  const auto threads = static_cast<unsigned>(args.num("threads", 0));
  const Outputs out{args.str("out", "bsar")};
  const std::string mode = args.str("mode", "geometric");
  if (mode != "geometric" && mode != "coherent" && mode != "hillshade") {
    throw std::invalid_argument("unknown --mode " + mode);
  }

  // ---------------------------------------------------------------- scene
  Scene scene;
  Heightmap dem;
  const std::string scene_spec = args.str("scene", "city");
  std::optional<EnuFrame> frame;
  if (args.has("origin")) {
    const auto o = args.nums("origin");
    if (o.size() != 3) {
      throw std::invalid_argument("--origin expects LAT,LON,H");
    }
    frame = EnuFrame({o[0], o[1], o[2]});
  }
  const bool has_dem = scene_spec.rfind("dem:", 0) == 0;
  Vec2d projected_offset{0.0, 0.0};
  std::filesystem::path dem_path;
  if (has_dem) {
    const std::filesystem::path path = scene_spec.substr(4);
    dem_path = path;
    if (path.extension() == ".asc") {
      dem = load_esri_ascii(path, args.has("geographic"));
    } else {
#ifdef BSAR_HAVE_GDAL
      dem = load_dem_gdal(path);
#else
      throw std::invalid_argument("only .asc DEMs are supported without BSAR_WITH_GDAL");
#endif
    }
    if (!dem.geographic) {
      // Projected grids (e.g. UTM) have coordinates ~1e6 m, where a float vertex
      // resolves only ~0.5 m. Re-centre the grid on the origin; footprints in
      // the same CRS get the same shift below.
      projected_offset = {dem.x0 + 0.5 * dem.dx * static_cast<double>(dem.nx - 1),
                          dem.y0 + 0.5 * dem.dy * static_cast<double>(dem.ny - 1)};
      dem.x0 -= projected_offset.x;
      dem.y0 -= projected_offset.y;
      std::printf("projected DEM re-centred: scene origin = (%.3f, %.3f) in the DEM's CRS\n",
                  projected_offset.x, projected_offset.y);
    }
    if (dem.geographic && !frame) {
      // Default origin: DEM centre at its mean height.
      double sum = 0.0;
      std::size_t n = 0;
      for (const float z : dem.z) {
        if (z != dem.nodata) {
          sum += z;
          ++n;
        }
      }
      frame = EnuFrame({dem.y0 + 0.5 * dem.dy * static_cast<double>(dem.ny - 1),
                        dem.x0 + 0.5 * dem.dx * static_cast<double>(dem.nx - 1),
                        n > 0 ? sum / static_cast<double>(n) : 0.0});
    }
    if (args.has("classes") != args.has("material-table")) {
      throw std::invalid_argument("--classes and --material-table go together");
    }
    if (args.has("classes")) {
      // Per-sample land-cover class -> material; unknown classes get dry soil.
      const Heightmap classes = load_esri_ascii(args.str("classes", ""), dem.geographic);
      if (classes.nx != dem.nx || classes.ny != dem.ny) {
        throw std::invalid_argument("--classes grid does not match the DEM grid");
      }
      const auto table = load_material_table(args.str("material-table", ""));
      const std::uint16_t fallback = scene.material_id("dry_soil");
      std::map<long, std::uint16_t> ids;
      std::vector<std::uint16_t> sample_material(classes.z.size(), fallback);
      for (std::size_t k = 0; k < classes.z.size(); ++k) {
        const long code = std::lround(classes.z[k]);
        auto it = ids.find(code);
        if (it == ids.end()) {
          const auto entry = table.find(static_cast<int>(code));
          it =
              ids.emplace(code, entry == table.end() ? fallback : scene.add_material(entry->second))
                  .first;
        }
        sample_material[k] = it->second;
      }
      scene.mesh = heightmap_to_mesh(dem, sample_material, frame ? &*frame : nullptr);
      std::printf("materials: %zu classes mapped to %zu materials\n", ids.size(),
                  scene.materials.size());
    } else {
      scene.mesh = heightmap_to_mesh(dem, scene.material_id("dry_soil"), frame ? &*frame : nullptr);
    }
  } else {
    scene = make_synthetic_scene(scene_spec);
  }
  if (args.has("buildings")) {
    FootprintSet fps = load_footprints(args.str("buildings", ""));
    if (!fps.geographic) {
      for (Footprint& fp : fps.footprints) {
        for (Vec2d& q : fp.ring) {
          q.x -= projected_offset.x;
          q.y -= projected_offset.y;
        }
      }
    }
    if (fps.geographic && !frame) {
      throw std::invalid_argument("geographic footprints need --origin");
    }
    extrude_footprints(
        fps, has_dem ? &dem : nullptr, frame ? &*frame : nullptr,
        [&](const std::string& m) { return scene.material_id(m); }, scene.mesh);
  }
  scene.frame = frame;
  scene.build();

  // ----------------------------------------------------------- trajectory
  std::unique_ptr<Trajectory> trajectory;
  const std::string platform = args.str("platform", "spaceborne");
  const bool airborne = platform == "airborne";
  if (!airborne && platform != "spaceborne") {
    throw std::invalid_argument("unknown --platform " + platform);
  }
  const double incidence = args.num("incidence", 35.0);
  const std::string look = args.str("look", "right");
  if (args.has("orbit")) {
    if (!frame) {
      throw std::invalid_argument("--orbit needs --origin (or a geographic DEM)");
    }
    Orbit orbit = load_eof(args.str("orbit", ""));
    if (args.has("time")) {
      orbit = orbit.window(parse_utc(args.str("time", "")), 600.0);
    }
    std::printf("orbit: %s .. %s UTC\n", format_utc(orbit.epoch()).c_str(),
                format_utc(orbit.epoch() + orbit.t_end()).c_str());
    trajectory = std::make_unique<OrbitTrajectory>(std::move(orbit), *frame);
  } else {
    trajectory = std::make_unique<LinearTrajectory>(make_side_looking_track(
        scene.bounds(), airborne ? 6000.0 : 693e3, airborne ? 150.0 : 7600.0, incidence,
        look == "left" ? LookSide::Left : LookSide::Right));
  }

  const Polarization pol = args.str("pol", "VV") == "HH" ? Polarization::HH : Polarization::VV;
  std::printf("scene: %zu triangles, bounds [%.0f %.0f %.0f] .. [%.0f %.0f %.0f] m\n",
              scene.mesh.num_triangles(), scene.bounds().lo.x, scene.bounds().lo.y,
              scene.bounds().lo.z, scene.bounds().hi.x, scene.bounds().hi.y, scene.bounds().hi.z);

  // -------------------------------------------------------------- render
  if (mode == "hillshade") {
    HillshadeConfig cfg;
    cfg.threads = threads;
    cfg.pixel_size = args.num("pixel", 0.0);
    cfg.sun_azimuth_deg = args.num("sun-azimuth", 315.0);
    cfg.sun_elevation_deg = args.num("sun-elevation", 45.0);
    args.reject_unknown();
    const Image<float> img = render_hillshade(scene, cfg);
    write_png(out.path("_hillshade.png"), to_u8(img, 0.0, 1.0));
    write_npy(out.path("_hillshade.npy"), img);
    std::printf("hillshade %zux%zu in %.2f s\n", img.width, img.height, seconds_since(t0));
    return 0;
  }

  if (mode == "geometric") {
    GeometricConfig cfg;
    cfg.range_spacing = args.num("range-spacing", 5.0);
    cfg.azimuth_spacing = args.num("azimuth-spacing", 5.0);
    cfg.rays_per_bin = args.num("rays-per-bin", 4.0);
    cfg.max_bounces = static_cast<int>(args.num("bounces", 3));
    cfg.polarization = pol;
    cfg.threads = threads;
    const bool do_geocode = args.has("geocode");
    if (do_geocode && !has_dem) {
      throw std::invalid_argument("--geocode needs a DEM scene (dem:PATH)");
    }
    args.reject_unknown();
    const auto t_render = std::chrono::steady_clock::now();
    const GeometricImage img = render_geometric(scene, *trajectory, cfg);
    const double render_s = seconds_since(t_render);

    write_npy(out.path("_intensity.npy"), img.intensity);
    write_png(out.path("_intensity.png"), to_db_u8_auto(img.intensity, 35.0));
    write_npy(out.path("_layover.npy"), img.layover);
    write_npy(out.path("_shadow.npy"), img.shadow);
    write_overlay(out, img);
    for (std::size_t b = 0; b < img.bounce.size(); ++b) {
      write_png(out.path("_bounce" + std::to_string(b + 1) + ".png"),
                to_db_u8_auto(img.bounce[b], 35.0));
    }
    std::ofstream meta(out.path("_meta.json"));
    meta << "{\n  \"mode\": \"geometric\",\n  \"range_bins\": " << img.intensity.width
         << ",\n  \"azimuth_lines\": " << img.intensity.height
         << ",\n  \"near_range_m\": " << img.near_range
         << ",\n  \"range_spacing_m\": " << img.range_spacing
         << ",\n  \"t_start_s\": " << img.t_start
         << ",\n  \"line_interval_s\": " << img.line_interval
         << ",\n  \"rays_traced\": " << img.rays_traced << ",\n  \"render_seconds\": " << render_s
         << "\n}\n";
    std::printf("geometric %zu lines x %zu bins, %zu rays in %.2f s (%.2f Mrays/s primary)\n",
                img.intensity.height, img.intensity.width, img.rays_traced, render_s,
                static_cast<double>(img.rays_traced) / render_s * 1e-6);

    if (do_geocode) {
      // Products on the DEM's own grid (same georeferencing as the input).
      const EnuFrame* enu = frame.has_value() ? &*frame : nullptr;
      if (dem.geographic && enu == nullptr) {
        throw std::logic_error("geocode: geographic DEM without an ENU frame");
      }
      std::vector<Vec3d> points(dem.nx * dem.ny);
      for (std::size_t j = 0; j < dem.ny; ++j) {
        for (std::size_t i = 0; i < dem.nx; ++i) {
          const double x = dem.x0 + static_cast<double>(i) * dem.dx;
          const double y = dem.y0 + static_cast<double>(j) * dem.dy;
          Vec3d& p = points[j * dem.nx + i];
          if (!dem.valid(i, j)) {
            p = {0.0, 0.0, std::numeric_limits<double>::quiet_NaN()};
          } else if (dem.geographic) {
            p = enu->geodetic_to_enu({y, x, dem.at(i, j)});
          } else {
            p = {x, y, dem.at(i, j)};
          }
        }
      }
      const auto t_geo = std::chrono::steady_clock::now();
      const GeocodedProducts geo =
          geocode(img, *trajectory, scene, points, dem.nx, dem.ny, threads);
      Heightmap grid = dem;  // georeferencing of the input file
      grid.x0 += projected_offset.x;
      grid.y0 += projected_offset.y;
      constexpr float kNoData = -9999.0f;
      auto as_float = [&](const auto& image, bool mask_invalid) {
        std::vector<float> v(image.data.size());
        for (std::size_t k = 0; k < v.size(); ++k) {
          v[k] =
              mask_invalid && geo.valid.data[k] == 0 ? kNoData : static_cast<float>(image.data[k]);
        }
        return v;
      };
      const std::filesystem::path prj = std::filesystem::path(dem_path).replace_extension(".prj");
      for (const auto& [name, values] : {std::pair{"beta0", as_float(geo.beta0, true)},
                                         std::pair{"layover", as_float(geo.layover, true)},
                                         std::pair{"shadow", as_float(geo.shadow, true)},
                                         std::pair{"valid", as_float(geo.valid, false)}}) {
        const std::filesystem::path file = out.path(std::string("_geo_") + name + ".asc");
        write_esri_ascii(file, grid, values, kNoData);
        if (std::filesystem::exists(prj)) {
          std::filesystem::copy_file(prj, std::filesystem::path(file).replace_extension(".prj"),
                                     std::filesystem::copy_options::overwrite_existing);
        }
      }
      std::printf("geocoded onto the DEM grid (%zu x %zu) in %.2f s\n", dem.nx, dem.ny,
                  seconds_since(t_geo));
    }
    return 0;
  }

  // coherent
  RadarParams radar;
  if (airborne) {
    radar.bandwidth = 100e6;
    radar.pulse_duration = 5e-6;
    radar.sample_rate = 150e6;
    radar.antenna_length = 2.0;
    radar.prf = 4.0 * 150.0 / radar.antenna_length;
  }
  ScattererConfig sc;
  sc.density = args.num("density", airborne ? 4.0 : 0.2);
  sc.seed = static_cast<std::uint64_t>(args.num("seed", 1));
  sc.polarization = pol;
  const double slant_res = slant_range_resolution(radar.bandwidth);
  const double pixel = args.num("pixel", 0.5 * std::min(slant_res / std::sin(deg_to_rad(incidence)),
                                                        radar.antenna_length / 2.0));
  args.reject_unknown();

  const EchoWindow window = auto_echo_window(scene.bounds(), *trajectory, radar);
  const Vec3d center(scene.bounds().center());
  const double t_c = find_zero_doppler_time(*trajectory, center);
  const FacetScatteringModel model;
  const auto scatterers = sample_scatterers(scene, trajectory->position(t_c), sc, model);
  std::printf("coherent: %zu scatterers, %.0f pulses\n", scatterers.size(),
              std::floor((window.t_end - window.t_start) * radar.prf) + 1);
  const RawData raw = synthesize_raw(scatterers, *trajectory, radar, window, threads);
  const CompressedData comp = range_compress(raw, 8, threads);

  // Ground grid over the scene, draped on the top surface.
  const Aabb& b = scene.bounds();
  PixelGrid grid;
  grid.width = static_cast<std::size_t>(std::ceil(b.extent().x / pixel));
  grid.height = static_cast<std::size_t>(std::ceil(b.extent().y / pixel));
  grid.origin = {b.lo.x + 0.5 * pixel, b.hi.y - 0.5 * pixel, 0.0};
  grid.step_x = {pixel, 0.0, 0.0};
  grid.step_y = {0.0, -pixel, 0.0};  // row 0 = north
  grid.heights.resize(grid.width * grid.height);
  for (std::size_t j = 0; j < grid.height; ++j) {
    for (std::size_t i = 0; i < grid.width; ++i) {
      Vec3d p = grid.position(i, j);
      Ray down{Vec3f(Vec3d{p.x, p.y, b.hi.z + 1.0}), {0.0f, 0.0f, -1.0f}};
      Hit hit;
      grid.heights[j * grid.width + i] =
          scene.intersect(down, hit) ? b.hi.z + 1.0 - hit.t : static_cast<double>(b.lo.z);
    }
  }
  BackprojectionConfig bp;
  bp.threads = threads;
  const auto slc = backproject(comp, grid, bp);
  Image<float> intensity(slc.width, slc.height);
  for (std::size_t i = 0; i < slc.size(); ++i) {
    intensity.data[i] = std::norm(slc.data[i]);
  }
  write_npy(out.path("_slc.npy"), slc);
  write_png(out.path("_intensity.png"), to_db_u8_auto(intensity, 30.0));
  std::printf("coherent %zux%zu px at %.2f m in %.2f s\n", slc.width, slc.height, pixel,
              seconds_since(t0));
  return 0;
} catch (const std::exception& e) {
  std::fprintf(stderr, "bsar_render: %s\n\n%s", e.what(), kUsage);
  return 1;
}
