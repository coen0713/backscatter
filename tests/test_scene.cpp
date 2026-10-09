#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <filesystem>
#include <fstream>

#include "backscatter/math/constants.hpp"
#include "backscatter/scene/buildings.hpp"
#include "backscatter/scene/heightmap.hpp"
#include "backscatter/scene/material.hpp"
#include "backscatter/scene/scene.hpp"

using namespace bsar;
using Catch::Approx;

namespace {

std::filesystem::path temp_file(const std::string& name, const std::string& contents) {
  const auto path = std::filesystem::temp_directory_path() / ("bsar_test_" + name);
  std::ofstream(path) << contents;
  return path;
}

}  // namespace

TEST_CASE("ESRI ASCII grid loading", "[scene][dem]") {
  const auto path = temp_file("grid.asc",
                              "ncols 3\nnrows 2\nxllcorner 100\nyllcorner 200\ncellsize 10\n"
                              "NODATA_value -9999\n1 2 3\n4 -9999 6\n");
  const Heightmap hm = load_esri_ascii(path, false);
  REQUIRE(hm.nx == 3);
  REQUIRE(hm.ny == 2);
  CHECK(hm.x0 == Approx(105.0));
  CHECK(hm.y0 == Approx(205.0));
  // First text row is the northern one (j = 1).
  CHECK(hm.at(0, 1) == 1.0f);
  CHECK(hm.at(2, 0) == 6.0f);
  CHECK_FALSE(hm.valid(1, 0));
  CHECK(std::isnan(hm.sample(110.0, 210.0)));  // touches no-data
  const TriangleMesh mesh = heightmap_to_mesh(hm, 0);
  CHECK(mesh.num_vertices() == 6);
  CHECK(mesh.num_triangles() == 0);  // both cells touch the no-data sample
  std::filesystem::remove(path);
}

TEST_CASE("Heightmap mesh is upward facing and sampled bilinearly", "[scene][dem]") {
  Heightmap hm;
  hm.nx = 4;
  hm.ny = 3;
  hm.dx = 2.0;
  hm.dy = 3.0;
  hm.z.resize(12);
  for (std::size_t j = 0; j < 3; ++j) {
    for (std::size_t i = 0; i < 4; ++i) {
      hm.at(i, j) = static_cast<float>(0.5 * static_cast<double>(i) * 2.0);  // z = 0.5 x
    }
  }
  CHECK(hm.sample(3.0, 4.0) == Approx(1.5));
  const TriangleMesh mesh = heightmap_to_mesh(hm, 0);
  CHECK(mesh.num_triangles() == 2 * 3 * 2);
  for (std::size_t t = 0; t < mesh.num_triangles(); ++t) {
    const Vec3f n = mesh.normal(t);
    CHECK(n.z > 0.0f);
    CHECK(n.x == Approx(-0.5 / std::sqrt(1.25)).margin(1e-6));
  }
}

TEST_CASE("Geographic DEMs are placed in the ENU frame", "[scene][dem]") {
  Heightmap hm;
  hm.nx = hm.ny = 2;
  hm.geographic = true;
  hm.x0 = 7.9;
  hm.y0 = 46.5;
  hm.dx = hm.dy = 0.001;
  hm.z = {1000.0f, 1000.0f, 1000.0f, 1000.0f};
  const EnuFrame frame({46.5, 7.9, 1000.0});
  const TriangleMesh mesh = heightmap_to_mesh(hm, 0, &frame);
  CHECK(length(mesh.vertex(0)) < 1e-3f);
  CHECK(mesh.vertex(3).y == Approx(111.2).margin(0.5));  // 0.001 deg of latitude
  CHECK_THROWS(heightmap_to_mesh(hm, 0, nullptr));
}

TEST_CASE("Polygon triangulation and footprint extrusion", "[scene][buildings]") {
  // L-shaped footprint, clockwise on purpose.
  const std::vector<Vec2d> ring{{0, 0}, {0, 20}, {10, 20}, {10, 10}, {20, 10}, {20, 0}};
  CHECK(signed_area(ring) == Approx(-300.0));
  const auto tris = triangulate_polygon(ring);
  REQUIRE(tris.size() == 3 * 4);
  double area = 0.0;
  for (std::size_t k = 0; k < tris.size(); k += 3) {
    const std::vector<Vec2d> tri{ring[tris[k]], ring[tris[k + 1]], ring[tris[k + 2]]};
    const double a = signed_area(tri);
    CHECK(a > 0.0);
    area += a;
  }
  CHECK(area == Approx(300.0));

  TriangleMesh mesh;
  extrude_footprint(ring, 5.0, 25.0, 0, mesh);
  CHECK(mesh.num_triangles() == 2 * 6 + 4);
  // Closed surface (open at the floor): outward walls, upward roof.
  const Aabb b = mesh.bounds();
  CHECK(b.lo.z == 5.0f);
  CHECK(b.hi.z == 25.0f);
  double roof = 0.0;
  for (std::size_t t = 0; t < mesh.num_triangles(); ++t) {
    const Vec3f n = mesh.normal(t);
    const auto [a, bb, c] = mesh.triangle(t);
    const Vec3f centroid = (a + bb + c) / 3.0f;
    if (n.z > 0.99f) {
      roof += mesh.area(t);
    } else {
      // Outward: stepping along the normal leaves the footprint.
      const Vec3f out = centroid + n * 0.5f;
      const bool inside_l = (out.x > 0 && out.x < 20 && out.y > 0 && out.y < 10) ||
                            (out.x > 0 && out.x < 10 && out.y > 0 && out.y < 20);
      CHECK_FALSE(inside_l);
    }
  }
  CHECK(roof == Approx(300.0));
}

TEST_CASE("Footprint file loading", "[scene][buildings]") {
  const auto path = temp_file("fp.txt",
                              "# test\ncrs enu\nbuilding 12.5 metal\n0 0\n10 0\n10 10\n0 10\nend\n"
                              "building 30\n50 50\n60 50\n55 60\nend\n");
  const FootprintSet set = load_footprints(path);
  CHECK_FALSE(set.geographic);
  REQUIRE(set.footprints.size() == 2);
  CHECK(set.footprints[0].height == 12.5);
  CHECK(set.footprints[0].material == "metal");
  CHECK(set.footprints[1].material == "concrete");
  CHECK(set.footprints[1].ring.size() == 3);

  Heightmap ground;
  ground.nx = ground.ny = 2;
  ground.x0 = ground.y0 = -100.0;
  ground.dx = ground.dy = 200.0;
  ground.z = {7.0f, 7.0f, 7.0f, 7.0f};
  Scene scene;
  extrude_footprints(
      set, &ground, nullptr, [&](const std::string& name) { return scene.material_id(name); },
      scene.mesh);
  CHECK(scene.mesh.bounds().lo.z == 7.0f);
  CHECK(scene.materials.size() == 2);
  std::filesystem::remove(path);
}

TEST_CASE("Fresnel coefficients", "[scene][material]") {
  Material m = materials::dry_soil();
  m.permittivity = {4.0, 0.0};
  // Normal incidence: (1 - sqrt(eps)) / (1 + sqrt(eps)) = -1/3 for both pols (H sign convention).
  CHECK(std::abs(fresnel_reflection(m, 1.0, Polarization::HH)) == Approx(1.0 / 3.0));
  CHECK(std::abs(fresnel_reflection(m, 1.0, Polarization::VV)) == Approx(1.0 / 3.0));
  // Brewster angle for eps = 4: tan(theta_B) = 2, |R_V| = 0.
  const double brewster = std::atan(2.0);
  CHECK(std::abs(fresnel_reflection(m, std::cos(brewster), Polarization::VV)) < 1e-9);
  // Grazing incidence reflects everything.
  CHECK(std::abs(fresnel_reflection(m, 0.0, Polarization::HH)) == Approx(1.0));
  CHECK(std::abs(fresnel_reflection(materials::metal(), 0.3, Polarization::HH)) == 1.0);
}

TEST_CASE("Facet scattering model", "[scene][material]") {
  const FacetScatteringModel model;
  const Material soil = materials::dry_soil();
  // Away from normal incidence, backscatter falls with incidence angle.
  double prev = model.backscatter(soil, std::cos(deg_to_rad(20.0)), Polarization::VV);
  for (double inc = 25.0; inc <= 70.0; inc += 5.0) {
    const double s = model.backscatter(soil, std::cos(deg_to_rad(inc)), Polarization::VV);
    CHECK(s < prev);
    prev = s;
  }
  const double s35 =
      10.0 * std::log10(model.backscatter(soil, std::cos(deg_to_rad(35.0)), Polarization::VV));
  CHECK(s35 > -16.0);
  CHECK(s35 < -8.0);
  // Calm water is dark off-nadir but bright at normal incidence.
  const Material water = materials::water();
  CHECK(model.backscatter(water, std::cos(deg_to_rad(35.0)), Polarization::VV) < 0.01);
  CHECK(model.backscatter(water, 1.0, Polarization::VV) > 10.0);
  // Bistatic reduces to monostatic.
  const Vec3d n{0, 0, 1};
  const Vec3d s = normalize(Vec3d{0.3, 0.1, 0.9});
  CHECK(model.bistatic(soil, n, s, s, Polarization::HH) ==
        Approx(model.backscatter(soil, dot(n, s), Polarization::HH)));
  CHECK(model.backscatter(soil, -0.1, Polarization::HH) == 0.0);
  CHECK_THROWS(materials::by_name("cheese"));
}

TEST_CASE("Synthetic scenes build", "[scene]") {
  for (const char* name : {"plane", "ridge", "dihedral", "city", "mountains"}) {
    const Scene scene = make_synthetic_scene(name);
    CHECK(scene.built());
    CHECK(scene.mesh.num_triangles() > 0);
    CHECK_FALSE(scene.bounds().empty());
  }
  CHECK_THROWS(make_synthetic_scene("atlantis"));
}
