// bsar_bench: reproducible performance numbers, printed as Markdown.
// Run through scripts/run_benchmarks.py, which records the commit and machine.

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "backscatter/eval/theory_checks.hpp"
#include "backscatter/geometry/bvh.hpp"
#include "backscatter/geometry/wide_bvh.hpp"
#include "backscatter/image/backprojection.hpp"
#include "backscatter/math/rng.hpp"
#include "backscatter/scene/heightmap.hpp"
#include "backscatter/scene/scene.hpp"
#include "backscatter/sensor/sar_geometry.hpp"
#include "backscatter/trace/coherent.hpp"
#include "backscatter/trace/geometric.hpp"
#include "backscatter/trace/sensor_rays.hpp"
#include "backscatter/util/cpu.hpp"
#include "backscatter/util/parallel.hpp"

#ifdef BSAR_HAVE_EMBREE
#include <embree4/rtcore.h>
#endif

using namespace bsar;

namespace {

#ifdef BSAR_HAVE_EMBREE
/// Intel Embree 4 over the same triangles, as a reference point. Default
/// scene flags; rtcIntersect1 per ray.
class EmbreeAccel {
 public:
  EmbreeAccel(const TriangleMesh& mesh, RTCBuildQuality quality) {
    device_ = rtcNewDevice(nullptr);
    scene_ = rtcNewScene(device_);
    rtcSetSceneBuildQuality(scene_, quality);
    RTCGeometry geom = rtcNewGeometry(device_, RTC_GEOMETRY_TYPE_TRIANGLE);
    rtcSetGeometryBuildQuality(geom, quality);
    auto* v = static_cast<float*>(rtcSetNewGeometryBuffer(geom, RTC_BUFFER_TYPE_VERTEX, 0,
                                                          RTC_FORMAT_FLOAT3, 3 * sizeof(float),
                                                          mesh.num_vertices()));
    for (std::size_t i = 0; i < mesh.num_vertices(); ++i) {
      const Vec3f q = mesh.vertex(static_cast<std::uint32_t>(i));
      v[3 * i] = q.x;
      v[3 * i + 1] = q.y;
      v[3 * i + 2] = q.z;
    }
    auto* idx = static_cast<unsigned*>(
        rtcSetNewGeometryBuffer(geom, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3,
                                3 * sizeof(unsigned), mesh.num_triangles()));
    for (std::size_t t = 0; t < mesh.num_triangles(); ++t) {
      const auto tri = mesh.indices(t);
      idx[3 * t] = tri[0];
      idx[3 * t + 1] = tri[1];
      idx[3 * t + 2] = tri[2];
    }
    rtcCommitGeometry(geom);
    rtcAttachGeometry(scene_, geom);
    rtcReleaseGeometry(geom);
    rtcCommitScene(scene_);
  }
  EmbreeAccel(const EmbreeAccel&) = delete;
  EmbreeAccel& operator=(const EmbreeAccel&) = delete;
  EmbreeAccel(EmbreeAccel&&) = delete;
  EmbreeAccel& operator=(EmbreeAccel&&) = delete;
  ~EmbreeAccel() {
    rtcReleaseScene(scene_);
    rtcReleaseDevice(device_);
  }

  bool intersect(const Ray& ray, Hit& hit) const {
    RTCRayHit rh{};
    rh.ray.org_x = ray.origin.x;
    rh.ray.org_y = ray.origin.y;
    rh.ray.org_z = ray.origin.z;
    rh.ray.dir_x = ray.dir.x;
    rh.ray.dir_y = ray.dir.y;
    rh.ray.dir_z = ray.dir.z;
    rh.ray.tnear = ray.tmin;
    rh.ray.tfar = ray.tmax;
    rh.ray.mask = 0xFFFFFFFFu;
    rh.hit.geomID = RTC_INVALID_GEOMETRY_ID;
    rtcIntersect1(scene_, &rh);
    if (rh.hit.geomID == RTC_INVALID_GEOMETRY_ID) {
      return false;
    }
    hit = {rh.ray.tfar, rh.hit.u, rh.hit.v, rh.hit.primID};
    return true;
  }

 private:
  RTCDevice device_ = nullptr;
  RTCScene scene_ = nullptr;
};
#endif

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
  std::printf("| Collapse SAH-16 to 4-wide | %.3f s | | %zu | |\n", collapse_s,
              wide.nodes().size());
  Bvh8 wide8;
  const double collapse8_s = time_best_of(quick ? 1 : 3, [&] { wide8.build(bvhs[0]); });
  std::printf("| Collapse SAH-16 to 8-wide | %.3f s | | %zu | |\n\n", collapse8_s,
              wide8.nodes().size());
  const char* wide8_kind = cpu_has_avx2_fma() ? "AVX2" : "scalar";

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
  std::printf("| 8-wide BVH (%s), SAH | 1 | %.2f |\n", wide8_kind, trace_rate(wide8, rays, 1));
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
  std::printf("| 4-wide BVH (SSE), SAH, coherent order | 1 | %.2f |\n", wide1_sorted);
  const double wide8_sorted = trace_rate(wide8, sorted, 1);
  std::printf("| 8-wide BVH (%s), SAH, coherent order | 1 | %.2f |\n", wide8_kind, wide8_sorted);
#ifdef BSAR_HAVE_EMBREE
  {
    std::unique_ptr<EmbreeAccel> embree;
    const double embree_build = time_best_of(
        1, [&] { embree = std::make_unique<EmbreeAccel>(scene.mesh, RTC_BUILD_QUALITY_MEDIUM); });
    const double e_random = trace_rate(*embree, rays, 1);
    const double e_sorted = trace_rate(*embree, sorted, 1);
    std::printf("| Embree %d.%d.%d (rtcIntersect1) | 1 | %.2f |\n", RTC_VERSION_MAJOR,
                RTC_VERSION_MINOR, RTC_VERSION_PATCH, e_random);
    std::printf("| Embree, coherent order | 1 | %.2f |\n\n", e_sorted);
    // Same closest hit? Compare hit distances on the first 100k rays.
    std::size_t agree = 0;
    std::size_t total = 0;
    for (std::size_t i = 0; i < std::min<std::size_t>(rays.size(), 100000); ++i) {
      Hit a;
      Hit b;
      const bool ha = wide8.intersect(rays[i], a);
      const bool hb = embree->intersect(rays[i], b);
      ++total;
      if (ha == hb && (!ha || std::abs(a.t - b.t) <= 1e-4f * std::max(1.0f, a.t))) {
        ++agree;
      }
    }
    const double best_sorted = std::max(wide1_sorted, wide8_sorted);
    const double best_random = std::max(wide1, trace_rate(wide8, rays, 1));
    std::printf(
        "Embree build: %.3f s (medium quality). Embree is %.2fx our best single-thread rate in "
        "coherent order (%.2fx in random order). Closest hits agree on %.3f%% of %zu rays.\n\n",
        embree_build, e_sorted / best_sorted, e_random / best_random,
        100.0 * static_cast<double>(agree) / static_cast<double>(total), total);
  }
#else
  std::printf("\n(Embree comparison not built: configure with -DBSAR_WITH_EMBREE=ON.)\n\n");
#endif
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
  blocked.kernel = BackprojectionKernel::Blocked;
  blocked.threads = 1;
  const double t_blocked =
      time_best_of(quick ? 1 : 3, [&] { (void)backproject(comp, pg, blocked); });
  std::printf("| Naive (pixel-major) | 1 | %.2f | 1.00x |\n", work / t_naive * 1e-6);
  std::printf("| Cache-blocked | 1 | %.2f | %.2fx |\n", work / t_blocked * 1e-6,
              t_naive / t_blocked);
  BackprojectionConfig best;
  best.kernel = resolve_kernel(BackprojectionKernel::Auto);
  const char* best_name =
      best.kernel == BackprojectionKernel::Simd ? "AVX2 SIMD, cache-blocked" : "Cache-blocked";
  if (best.kernel == BackprojectionKernel::Simd) {
    best.threads = 1;
    const double t_simd = time_best_of(quick ? 1 : 3, [&] { (void)backproject(comp, pg, best); });
    std::printf("| %s | 1 | %.2f | %.2fx |\n", best_name, work / t_simd * 1e-6, t_naive / t_simd);
  }
  best.threads = hw;
  const double t_mt = time_best_of(quick ? 1 : 3, [&] { (void)backproject(comp, pg, best); });
  std::printf("| %s | %u | %.2f | %.2fx |\n", best_name, hw, work / t_mt * 1e-6, t_naive / t_mt);
  return 0;
}
