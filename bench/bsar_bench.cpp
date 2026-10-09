// bsar_bench: reproducible performance numbers, printed as Markdown.
// Run through scripts/run_benchmarks.py, which records the commit and machine.

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <functional>
#include <string>
#include <vector>

#include "backscatter/eval/theory_checks.hpp"
#include "backscatter/geometry/bvh.hpp"
#include "backscatter/geometry/bvh4.hpp"
#include "backscatter/image/backprojection.hpp"
#include "backscatter/math/rng.hpp"
#include "backscatter/scene/heightmap.hpp"
#include "backscatter/scene/scene.hpp"
#include "backscatter/sensor/sar_geometry.hpp"
#include "backscatter/trace/coherent.hpp"
#include "backscatter/trace/geometric.hpp"
#include "backscatter/trace/sensor_rays.hpp"
#include "backscatter/util/parallel.hpp"

using namespace bsar;

namespace {

double time_best_of(int repeats, const std::function<void()>& fn) {
  double best = 1e300;
  for (int r = 0; r < repeats; ++r) {
    const auto t0 = std::chrono::steady_clock::now();
    fn();
    best = std::min(best,
                    std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
  }
  return best;
}

/// Rays from a spaceborne sensor towards random terrain points, clipped to the
/// scene bounds: the same distribution the geometric integrator produces.
std::vector<Ray> make_sar_rays(const Scene& scene, std::size_t n) {
  const auto track = make_side_looking_track(scene.bounds(), 693e3, 7600.0, 35.0, LookSide::Right);
  const Aabb& b = scene.bounds();
  const CounterRng rng(17);
  std::vector<Ray> rays;
  rays.reserve(n);
  for (std::uint64_t i = 0; rays.size() < n; ++i) {
    const auto u = rng.uniform4(i, 0);
    const Vec3d target{b.lo.x + u[0] * (b.hi.x - b.lo.x), b.lo.y + u[1] * (b.hi.y - b.lo.y),
                       b.lo.z + u[2] * (b.hi.z - b.lo.z)};
    const Vec3d p = track.position(target.y / 7600.0);
    if (auto c = clip_to_bounds(p, normalize(target - p), b, 1.0)) {
      rays.push_back(c->ray);
    }
  }
  return rays;
}

template <typename Accel>
double trace_rate(const Accel& accel, const std::vector<Ray>& rays, unsigned threads,
                  std::size_t* hits_out = nullptr) {
  std::vector<unsigned char> hit(rays.size());
  const double s = time_best_of(3, [&] {
    parallel_for(
        0, rays.size(), threads,
        [&](std::size_t i) {
          Hit h;
          hit[i] = accel.intersect(rays[i], h) ? 1 : 0;
        },
        4096);
  });
  if (hits_out != nullptr) {
    *hits_out = 0;
    for (const auto x : hit) {
      *hits_out += x;
    }
  }
  return static_cast<double>(rays.size()) / s * 1e-6;
}

}  // namespace

int main(int argc, char** argv) {
  const bool quick = argc > 1 && std::strcmp(argv[1], "--quick") == 0;
  const std::size_t grid = quick ? 257 : 1025;
  const std::size_t n_rays = quick ? 200'000 : 2'000'000;
  const unsigned hw = resolve_threads(0);

  // ------------------------------------------------------------- scene
  const Heightmap hm = make_fractal_terrain(grid, 30.0, 2000.0, 7);
  Scene scene;
  scene.mesh = heightmap_to_mesh(hm, scene.material_id("dry_soil"));
  scene.build();
  const std::size_t tris = scene.mesh.num_triangles();
  std::printf("## BVH construction\n\n");
  std::printf("Fractal terrain, %zu x %zu samples at 30 m (%zu triangles).\n\n", grid, grid, tris);
  std::printf("| Builder | Build time | SAH cost | Nodes | Max depth |\n|---|---|---|---|---|\n");

  struct Variant {
    const char* name;
    BvhBuildConfig cfg;
  };
  BvhBuildConfig sah16;
  BvhBuildConfig sah32;
  sah32.num_bins = 32;
  BvhBuildConfig median;
  median.strategy = BvhBuildStrategy::ObjectMedian;
  const std::vector<Variant> variants{
      {"Binned SAH, 16 bins", sah16}, {"Binned SAH, 32 bins", sah32}, {"Object median", median}};
  std::vector<Bvh> bvhs(variants.size());
  for (std::size_t v = 0; v < variants.size(); ++v) {
    const double s =
        time_best_of(quick ? 1 : 3, [&] { bvhs[v].build(scene.mesh, variants[v].cfg); });
    const BvhStats& st = bvhs[v].stats();
    std::printf("| %s | %.3f s | %.2f | %zu | %zu |\n", variants[v].name, s, st.sah_cost,
                st.num_nodes, st.max_depth);
  }
  Bvh4 wide;
  const double collapse_s = time_best_of(quick ? 1 : 3, [&] { wide.build(bvhs[0]); });
  std::printf("| Collapse SAH-16 to 4-wide | %.3f s | | %zu | |\n\n", collapse_s,
              wide.nodes().size());

  // ---------------------------------------------------------- traversal
  const std::vector<Ray> rays = make_sar_rays(scene, n_rays);
  std::size_t hits = 0;
  std::printf("## Ray traversal (closest hit)\n\n");
  std::printf("%zu SAR-geometry rays (35 deg incidence from 693 km, clipped to the scene).\n\n",
              rays.size());
  std::printf("| Structure | Threads | Mrays/s |\n|---|---|---|\n");
  std::printf("| Binary BVH, object median | 1 | %.2f |\n", trace_rate(bvhs[2], rays, 1));
  std::printf("| Binary BVH, SAH | 1 | %.2f |\n", trace_rate(bvhs[0], rays, 1));
  const double wide1 = trace_rate(wide, rays, 1, &hits);
  std::printf("| 4-wide BVH (SSE), SAH | 1 | %.2f |\n", wide1);
  // The integrator issues rays line by line, sweeping the look angle, so
  // neighbouring rays touch neighbouring nodes. Same rays in that order:
  std::vector<Ray> sorted = rays;
  std::sort(sorted.begin(), sorted.end(), [](const Ray& a, const Ray& b) {
    const auto la = static_cast<long>(std::floor(a.origin.y / 30.0f));
    const auto lb = static_cast<long>(std::floor(b.origin.y / 30.0f));
    return la != lb ? la < lb : a.origin.x < b.origin.x;
  });
  std::printf("| Binary BVH, SAH, coherent order | 1 | %.2f |\n", trace_rate(bvhs[0], sorted, 1));
  const double wide1_sorted = trace_rate(wide, sorted, 1);
  std::printf("| 4-wide BVH (SSE), SAH, coherent order | 1 | %.2f |\n\n", wide1_sorted);
  std::printf(
      "Hit rate: %.1f%%. Random order is the cache-hostile worst case; the renderer "
      "traces in coherent order.\n\n",
      100.0 * static_cast<double>(hits) / static_cast<double>(rays.size()));

  std::printf("## Thread scaling (4-wide BVH)\n\n");
  std::printf("| Threads | Mrays/s | Speed-up | Efficiency |\n|---|---|---|---|\n");
  for (unsigned t = 1; t <= hw; t *= 2) {
    const double rate = t == 1 ? wide1 : trace_rate(wide, rays, t);
    std::printf("| %u | %.2f | %.2fx | %.0f%% |\n", t, rate, rate / wide1,
                100.0 * rate / wide1 / t);
    if (t * 2 > hw && t != hw) {
      const double r = trace_rate(wide, rays, hw);
      std::printf("| %u | %.2f | %.2fx | %.0f%% |\n", hw, r, r / wide1, 100.0 * r / wide1 / hw);
    }
  }
  std::printf("\n");

  // --------------------------------------------------- geometric render
  const auto track = make_side_looking_track(scene.bounds(), 693e3, 7600.0, 35.0, LookSide::Right);
  GeometricConfig gcfg;
  gcfg.range_spacing = 10.0;
  gcfg.azimuth_spacing = quick ? 60.0 : 15.0;
  gcfg.max_bounces = 3;
  GeometricImage gimg;
  const double gs = time_best_of(1, [&] { gimg = render_geometric(scene, track, gcfg); });
  std::printf("## Geometric render\n\n");
  std::printf(
      "Same terrain, 3 bounces, %zu lines x %zu bins, all %u threads: %zu primary rays "
      "in %.2f s (%.2f M primary rays/s, including shadow and bounce rays).\n\n",
      gimg.intensity.height, gimg.intensity.width, hw, gimg.rays_traced, gs,
      static_cast<double>(gimg.rays_traced) / gs * 1e-6);

  // ------------------------------------------------------- backprojection
  const RadarParams radar = theory_radar();
  const LinearTrajectory bp_track({-4000.0, 0.0, 3000.0}, {0.0, 100.0, 0.0}, -100.0, 100.0);
  const std::vector<Scatterer> pts{{{0.0, 0.0, 0.0}, 1.0}, {{5.0, 3.0, 0.0}, 0.7}};
  const EchoWindow w{-1.0, 1.0, 4900.0, 5100.0};
  const CompressedData comp = range_compress(synthesize_raw(pts, bp_track, radar, w), 8);
  PixelGrid pg;
  pg.width = pg.height = quick ? 64 : 256;
  pg.origin = {-16.0, -16.0, 0.0};
  pg.step_x = {32.0 / static_cast<double>(pg.width), 0.0, 0.0};
  pg.step_y = {0.0, 32.0 / static_cast<double>(pg.height), 0.0};
  const double work = static_cast<double>(pg.width * pg.height * comp.num_pulses);
  std::printf("## Backprojection\n\n");
  std::printf("%zu x %zu pixels x %zu pulses (%.1f M pixel-pulses).\n\n", pg.width, pg.height,
              comp.num_pulses, work * 1e-6);
  std::printf("| Kernel | Threads | M pixel-pulses/s | Speed-up |\n|---|---|---|---|\n");
  BackprojectionConfig naive;
  naive.kernel = BackprojectionKernel::Naive;
  naive.threads = 1;
  const double t_naive = time_best_of(quick ? 1 : 3, [&] { (void)backproject(comp, pg, naive); });
  BackprojectionConfig blocked;
  blocked.threads = 1;
  const double t_blocked =
      time_best_of(quick ? 1 : 3, [&] { (void)backproject(comp, pg, blocked); });
  blocked.threads = hw;
  const double t_mt = time_best_of(quick ? 1 : 3, [&] { (void)backproject(comp, pg, blocked); });
  std::printf("| Naive (pixel-major) | 1 | %.2f | 1.00x |\n", work / t_naive * 1e-6);
  std::printf("| Cache-blocked | 1 | %.2f | %.2fx |\n", work / t_blocked * 1e-6,
              t_naive / t_blocked);
  std::printf("| Cache-blocked | %u | %.2f | %.2fx |\n", hw, work / t_mt * 1e-6, t_naive / t_mt);
  return 0;
}
