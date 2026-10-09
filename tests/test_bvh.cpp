#include <algorithm>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <limits>

#include "backscatter/geometry/bvh.hpp"
#include "backscatter/geometry/triangle.hpp"
#include "backscatter/geometry/wide_bvh.hpp"
#include "backscatter/math/rng.hpp"
#include "backscatter/scene/heightmap.hpp"

using namespace bsar;
using Catch::Approx;

namespace {

TriangleMesh random_triangle_soup(std::size_t n, std::uint64_t seed) {
  const CounterRng rng(seed);
  TriangleMesh mesh;
  for (std::size_t t = 0; t < n; ++t) {
    const auto c = rng.uniform4(t, 0);
    const Vec3d center{c[0] * 100.0, c[1] * 100.0, c[2] * 100.0};
    const double size = 0.5 + 4.0 * c[3];
    std::uint32_t idx[3];
    for (std::uint64_t v = 0; v < 3; ++v) {
      const auto o = rng.uniform4(t, v + 1);
      idx[v] = mesh.add_vertex(Vec3f(center + Vec3d{o[0] - 0.5, o[1] - 0.5, o[2] - 0.5} * size));
    }
    mesh.add_triangle(idx[0], idx[1], idx[2], 0);
  }
  return mesh;
}

Ray random_ray(const CounterRng& rng, std::uint64_t i) {
  const auto a = rng.uniform4(i, 100);
  const auto b = rng.uniform4(i, 101);
  const Vec3d origin{a[0] * 140.0 - 20.0, a[1] * 140.0 - 20.0, a[2] * 140.0 - 20.0};
  const Vec3d target{b[0] * 100.0, b[1] * 100.0, b[2] * 100.0};
  Ray r;
  r.origin = Vec3f(origin);
  r.dir = Vec3f(normalize(target - origin));
  r.tmax = a[3] < 0.2 ? 30.0f : std::numeric_limits<float>::infinity();
  return r;
}

bool brute_force(const TriangleMesh& mesh, const Ray& ray, Hit& hit) {
  const WatertightRay w(ray);
  bool found = false;
  float tmax = ray.tmax;
  for (std::size_t t = 0; t < mesh.num_triangles(); ++t) {
    const auto [a, b, c] = mesh.triangle(t);
    float th = 0.0f;
    float u = 0.0f;
    float v = 0.0f;
    if (intersect_triangle(w, a, b, c, ray.tmin, tmax, th, u, v)) {
      tmax = th;
      hit = {th, u, v, static_cast<std::uint32_t>(t)};
      found = true;
    }
  }
  return found;
}

}  // namespace

TEST_CASE("Watertight intersection reports correct t and barycentrics", "[bvh][triangle]") {
  const Vec3f v0{0, 0, 0};
  const Vec3f v1{1, 0, 0};
  const Vec3f v2{0, 1, 0};
  Ray ray{{0.25f, 0.5f, 2.0f}, {0, 0, -1}};
  float t = 0.0f;
  float u = 0.0f;
  float v = 0.0f;
  REQUIRE(intersect_triangle(WatertightRay(ray), v0, v1, v2, 0.0f, 10.0f, t, u, v));
  CHECK(t == Approx(2.0f));
  CHECK(u == Approx(0.25f));
  CHECK(v == Approx(0.5f));
  ray.origin = {0.75f, 0.75f, 2.0f};
  CHECK_FALSE(intersect_triangle(WatertightRay(ray), v0, v1, v2, 0.0f, 10.0f, t, u, v));
}

TEST_CASE("Rays through shared edges and vertices never leak", "[bvh][triangle][watertight]") {
  // A terrain-like grid, hit at grazing angles exactly through vertices and
  // edge midpoints: every ray must hit something.
  Heightmap hm;
  hm.nx = hm.ny = 9;
  hm.dx = hm.dy = 1.0;
  hm.z.assign(81, 0.0f);
  for (std::size_t j = 0; j < 9; ++j) {
    for (std::size_t i = 0; i < 9; ++i) {
      hm.at(i, j) = static_cast<float>(0.1 * std::sin(static_cast<double>(i * 3 + j)));
    }
  }
  const TriangleMesh mesh = heightmap_to_mesh(hm, 0);
  Bvh bvh;
  bvh.build(mesh);
  Bvh4 bvh4;
  bvh4.build(bvh);
  Bvh8 bvh8;
  bvh8.build(bvh);
  std::size_t misses = 0;
  std::size_t total = 0;
  for (std::size_t t = 0; t < mesh.num_triangles(); ++t) {
    // Aim at every vertex and edge midpoint of interior triangles: the
    // points most likely to fall between neighbouring triangles.
    const auto [a, b, c] = mesh.triangle(t);
    const Vec3d centroid = (Vec3d(a) + Vec3d(b) + Vec3d(c)) / 3.0;
    if (centroid.x < 1.0 || centroid.x > 7.0 || centroid.y < 1.0 || centroid.y > 7.0) {
      continue;
    }
    // A segment from above the heightfield (z > 0.1) to below it (z < -0.1)
    // must cross the surface, so any miss is a leak. (A grazing ray that only
    // touches a convex vertex may legitimately miss, so the segment always
    // crosses the full height range.)
    for (const Vec3d target : {Vec3d(a), Vec3d(b), Vec3d(c), (Vec3d(a) + Vec3d(b)) * 0.5,
                               (Vec3d(b) + Vec3d(c)) * 0.5, (Vec3d(a) + Vec3d(c)) * 0.5}) {
      for (const double grazing_deg : {5.0, 10.0, 30.0, 60.0, 89.0}) {
        for (const double heading_deg : {0.0, 37.0, 90.0, 225.0}) {
          const double g = grazing_deg * 3.14159265358979 / 180.0;
          const double h = heading_deg * 3.14159265358979 / 180.0;
          const Vec3d dir{std::cos(g) * std::cos(h), std::cos(g) * std::sin(h), -std::sin(g)};
          const double half = 0.21 / std::sin(g);  // clears the +/-0.1 m relief
          const Vec3d start = target - dir * half;
          const Vec3d end = target + dir * half;
          if (std::min({start.x, start.y, end.x, end.y}) < 0.0 ||
              std::max({start.x, start.y, end.x, end.y}) > 8.0) {
            continue;  // segment would leave the heightfield's footprint
          }
          Ray ray;
          ray.origin = Vec3f(start);
          ray.dir = Vec3f(dir);
          ray.tmax = static_cast<float>(2.0 * half);
          Hit h1;
          Hit h2;
          ++total;
          if (!bvh.intersect(ray, h1) || !bvh4.intersect(ray, h2) || !bvh8.occluded(ray)) {
            ++misses;
          }
        }
      }
    }
  }
  CHECK(total > 0);
  CHECK(misses == 0);
}

TEST_CASE("BVH traversal matches brute force", "[bvh]") {
  const TriangleMesh mesh = random_triangle_soup(3000, 11);
  for (const auto strategy : {BvhBuildStrategy::BinnedSah, BvhBuildStrategy::ObjectMedian}) {
    BvhBuildConfig cfg;
    cfg.strategy = strategy;
    Bvh bvh;
    bvh.build(mesh, cfg);
    Bvh4 bvh4;
    bvh4.build(bvh);
    Bvh8 bvh8;
    bvh8.build(bvh);
    const CounterRng rng(5);
    std::size_t hits = 0;
    for (std::uint64_t i = 0; i < 4000; ++i) {
      const Ray ray = random_ray(rng, i);
      Hit ref;
      Hit a;
      Hit b;
      const bool expect = brute_force(mesh, ray, ref);
      REQUIRE(bvh.intersect(ray, a) == expect);
      REQUIRE(bvh4.intersect(ray, b) == expect);
      Hit c;
      REQUIRE(bvh8.intersect(ray, c) == expect);
      REQUIRE(bvh8.occluded(ray) == expect);
      REQUIRE(bvh.occluded(ray) == expect);
      REQUIRE(bvh4.occluded(ray) == expect);
      if (expect) {
        ++hits;
        REQUIRE(a.t == ref.t);
        REQUIRE(b.t == ref.t);
        REQUIRE(c.t == ref.t);
      }
    }
    CHECK(hits > 500);
  }
}

TEST_CASE("SAH build beats object median on SAH cost and respects leaf size", "[bvh]") {
  const TriangleMesh mesh = random_triangle_soup(20000, 3);
  Bvh sah;
  sah.build(mesh);
  BvhBuildConfig median_cfg;
  median_cfg.strategy = BvhBuildStrategy::ObjectMedian;
  Bvh median;
  median.build(mesh, median_cfg);
  CHECK(sah.stats().sah_cost < median.stats().sah_cost);
  std::size_t prims = 0;
  for (const BvhNode& n : sah.nodes()) {
    if (n.is_leaf()) {
      CHECK(n.count <= 4);
      prims += n.count;
    }
  }
  CHECK(prims == mesh.num_triangles());
  CHECK(sah.stats().max_depth < 64);
}

TEST_CASE("Degenerate inputs build valid hierarchies", "[bvh]") {
  TriangleMesh mesh;
  for (int i = 0; i < 100; ++i) {  // 100 identical triangles: no centroid spread
    const auto a = mesh.add_vertex({0, 0, 0});
    const auto b = mesh.add_vertex({1, 0, 0});
    const auto c = mesh.add_vertex({0, 1, 0});
    mesh.add_triangle(a, b, c, 0);
  }
  Bvh bvh;
  bvh.build(mesh);
  Bvh4 bvh4;
  bvh4.build(bvh);
  Ray ray{{0.2f, 0.2f, 1.0f}, {0, 0, -1}};
  Hit hit;
  CHECK(bvh.intersect(ray, hit));
  CHECK(bvh4.intersect(ray, hit));
  CHECK(hit.t == Approx(1.0f));

  Bvh empty;
  empty.build(TriangleMesh{});
  CHECK_FALSE(empty.intersect(ray, hit));
  CHECK_FALSE(empty.occluded(ray));
}
