#include "backscatter/scene/scene.hpp"

#include <cmath>
#include <stdexcept>

#include "backscatter/math/constants.hpp"
#include "backscatter/scene/buildings.hpp"
#include "backscatter/scene/heightmap.hpp"

namespace bsar {

std::uint16_t Scene::add_material(const Material& m) {
  for (std::size_t i = 0; i < materials.size(); ++i) {
    if (materials[i].name == m.name) {
      return static_cast<std::uint16_t>(i);
    }
  }
  materials.push_back(m);
  return static_cast<std::uint16_t>(materials.size() - 1);
}

std::uint16_t Scene::material_id(const std::string& name) {
  for (std::size_t i = 0; i < materials.size(); ++i) {
    if (materials[i].name == name) {
      return static_cast<std::uint16_t>(i);
    }
  }
  return add_material(materials::by_name(name));
}

void Scene::build(const BvhBuildConfig& config, AccelKind accel) {
  for (std::size_t t = 0; t < mesh.num_triangles(); ++t) {
    if (mesh.material(t) >= materials.size()) {
      throw std::runtime_error("Scene::build: triangle references a missing material");
    }
  }
  bvh_.build(mesh, config);
  accel_ = accel;
  if (accel_ == AccelKind::Wide4) {
    bvh4_.build(bvh_);
  }
  bounds_ = bvh_.bounds();
  built_ = true;
}

bool Scene::intersect(const Ray& ray, Hit& hit) const {
  return accel_ == AccelKind::Wide4 ? bvh4_.intersect(ray, hit) : bvh_.intersect(ray, hit);
}

bool Scene::occluded(const Ray& ray) const {
  return accel_ == AccelKind::Wide4 ? bvh4_.occluded(ray) : bvh_.occluded(ray);
}

namespace {

/// Extrude a z(x) profile along y over [-half_length, half_length].
void add_profile_strip(TriangleMesh& mesh, std::span<const Vec2d> profile_xz, double half_length,
                       std::uint16_t material) {
  const auto first = static_cast<std::uint32_t>(mesh.num_vertices());
  for (const Vec2d& p : profile_xz) {
    mesh.add_vertex(Vec3f(Vec3d{p.x, -half_length, p.y}));
    mesh.add_vertex(Vec3f(Vec3d{p.x, half_length, p.y}));
  }
  for (std::size_t i = 0; i + 1 < profile_xz.size(); ++i) {
    const auto a = first + static_cast<std::uint32_t>(2 * i);  // (x_i, south)
    const auto b = a + 1;                                      // (x_i, north)
    const auto c = a + 2;                                      // (x_i+1, south)
    const auto d = a + 3;                                      // (x_i+1, north)
    mesh.add_triangle(a, c, d, material);
    mesh.add_triangle(a, d, b, material);
  }
}

}  // namespace

Scene make_plane_scene(double size, const std::string& material) {
  Scene scene;
  const auto m = scene.material_id(material);
  const double h = 0.5 * size;
  const std::vector<Vec2d> profile{{-h, 0.0}, {h, 0.0}};
  add_profile_strip(scene.mesh, profile, h, m);
  scene.build();
  return scene;
}

Scene make_ridge_scene(double west_slope_deg, double east_slope_deg, double height,
                       double ground_size) {
  Scene scene;
  const auto m = scene.material_id("dry_soil");
  const double h = 0.5 * ground_size;
  const double west_run = height / std::tan(deg_to_rad(west_slope_deg));
  const double east_run = height / std::tan(deg_to_rad(east_slope_deg));
  if (west_run >= h || east_run >= h) {
    throw std::invalid_argument("make_ridge_scene: ridge does not fit on the ground");
  }
  const std::vector<Vec2d> profile{
      {-h, 0.0}, {-west_run, 0.0}, {0.0, height}, {east_run, 0.0}, {h, 0.0}};
  add_profile_strip(scene.mesh, profile, h, m);
  scene.build();
  return scene;
}

Scene make_dihedral_scene(double ground_size, double width, double depth, double height) {
  Scene scene;
  const auto ground = scene.material_id("dry_soil");
  const auto wall = scene.material_id("concrete");
  const double h = 0.5 * ground_size;
  const std::vector<Vec2d> profile{{-h, 0.0}, {h, 0.0}};
  add_profile_strip(scene.mesh, profile, h, ground);
  const double hx = 0.5 * width;
  const double hy = 0.5 * depth;
  const std::vector<Vec2d> ring{{-hx, -hy}, {hx, -hy}, {hx, hy}, {-hx, hy}};
  extrude_footprint(ring, -0.5, height, wall, scene.mesh);  // base slightly below ground
  scene.build();
  return scene;
}

Scene make_synthetic_scene(std::string_view name) {
  if (name == "plane") {
    return make_plane_scene(1000.0);
  }
  if (name == "ridge") {
    return make_ridge_scene(50.0, 30.0, 400.0, 3000.0);
  }
  if (name == "dihedral") {
    return make_dihedral_scene(400.0, 40.0, 40.0, 30.0);
  }
  if (name == "city") {
    Scene scene;
    const auto ground = scene.material_id("dry_soil");
    const auto concrete = scene.material_id("concrete");
    const auto metal = scene.material_id("metal");
    const std::vector<Vec2d> profile{{-400.0, 0.0}, {400.0, 0.0}};
    add_profile_strip(scene.mesh, profile, 400.0, ground);
    // A small grid of blocks of varying height, plus one rotated tower.
    for (int bx = -2; bx <= 2; ++bx) {
      for (int by = -2; by <= 2; ++by) {
        const double cx = bx * 120.0;
        const double cy = by * 120.0;
        const double hgt = 15.0 + 10.0 * (((bx + 2) * 5 + (by + 2)) % 6);
        const std::vector<Vec2d> ring{
            {cx - 30, cy - 25}, {cx + 30, cy - 25}, {cx + 30, cy + 25}, {cx - 30, cy + 25}};
        extrude_footprint(ring, -0.5, hgt, (bx + by) % 3 == 0 ? metal : concrete, scene.mesh);
      }
    }
    scene.build();
    return scene;
  }
  if (name == "mountains") {
    Scene scene;
    const Heightmap hm = make_fractal_terrain(257, 30.0, 1500.0, 7);
    scene.mesh = heightmap_to_mesh(hm, scene.material_id("dry_soil"));
    scene.build();
    return scene;
  }
  throw std::invalid_argument("unknown synthetic scene '" + std::string(name) +
                              "' (expected plane, ridge, dihedral, city or mountains)");
}

}  // namespace bsar
