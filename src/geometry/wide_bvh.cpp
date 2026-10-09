#include "backscatter/geometry/wide_bvh.hpp"

#include <algorithm>
#include <array>
#include <limits>

#include "backscatter/geometry/ray_box.hpp"
#include "backscatter/geometry/triangle.hpp"
#include "backscatter/util/cpu.hpp"

#if BSAR_X86_64
#include <immintrin.h>
#endif

namespace bsar {
namespace {

constexpr std::size_t kStackSize = 512;

template <int W>
void set_slot(WideBvhNode<W>& node, int slot, const Aabb& box, std::uint32_t child,
              std::uint32_t count) {
  const auto s = static_cast<std::size_t>(slot);
  node.lo_x[s] = box.lo.x;
  node.lo_y[s] = box.lo.y;
  node.lo_z[s] = box.lo.z;
  node.hi_x[s] = box.hi.x;
  node.hi_y[s] = box.hi.y;
  node.hi_z[s] = box.hi.z;
  node.child[s] = child;
  node.count[s] = count;
}

template <int W>
WideBvhNode<W> empty_node() {
  WideBvhNode<W> node{};
  for (int s = 0; s < W; ++s) {
    set_slot(node, s, Aabb{}, WideBvhNode<W>::kEmpty, 0);
  }
  return node;
}

template <int W>
unsigned empty_slots(const WideBvhNode<W>& node) {
  unsigned mask = 0;
  for (int s = 0; s < W; ++s) {
    if (node.child[static_cast<std::size_t>(s)] == WideBvhNode<W>::kEmpty) {
      mask |= 1u << s;
    }
  }
  return mask;
}

// ---------------------------------------------------------------- box tests
// Each returns a bit mask of the children whose boxes overlap [tmin, tmax]
// and writes their entry distances.

struct ScalarBoxes {
  template <int W>
  unsigned operator()(const WideBvhNode<W>& node, const Vec3f& org, const Vec3f& inv, float tmin,
                      float tmax, float* t_near) const {
    unsigned mask = 0;
    for (int s = 0; s < W; ++s) {
      const auto i = static_cast<std::size_t>(s);
      Aabb box;
      box.lo = {node.lo_x[i], node.lo_y[i], node.lo_z[i]};
      box.hi = {node.hi_x[i], node.hi_y[i], node.hi_z[i]};
      if (intersect_aabb(box, org, inv, tmin, tmax, t_near[i])) {
        mask |= 1u << s;
      }
    }
    return mask & ~empty_slots(node);
  }
};

#if BSAR_X86_64
struct SseBoxes {
  unsigned operator()(const WideBvhNode<4>& node, const Vec3f& org, const Vec3f& inv, float tmin,
                      float tmax, float* t_near) const {
    const __m128 ox = _mm_set1_ps(org.x);
    const __m128 oy = _mm_set1_ps(org.y);
    const __m128 oz = _mm_set1_ps(org.z);
    const __m128 ix = _mm_set1_ps(inv.x);
    const __m128 iy = _mm_set1_ps(inv.y);
    const __m128 iz = _mm_set1_ps(inv.z);
    const __m128 tx0 = _mm_mul_ps(_mm_sub_ps(_mm_load_ps(node.lo_x), ox), ix);
    const __m128 tx1 = _mm_mul_ps(_mm_sub_ps(_mm_load_ps(node.hi_x), ox), ix);
    const __m128 ty0 = _mm_mul_ps(_mm_sub_ps(_mm_load_ps(node.lo_y), oy), iy);
    const __m128 ty1 = _mm_mul_ps(_mm_sub_ps(_mm_load_ps(node.hi_y), oy), iy);
    const __m128 tz0 = _mm_mul_ps(_mm_sub_ps(_mm_load_ps(node.lo_z), oz), iz);
    const __m128 tz1 = _mm_mul_ps(_mm_sub_ps(_mm_load_ps(node.hi_z), oz), iz);
    const __m128 t0 = _mm_max_ps(_mm_max_ps(_mm_min_ps(tx0, tx1), _mm_min_ps(ty0, ty1)),
                                 _mm_max_ps(_mm_min_ps(tz0, tz1), _mm_set1_ps(tmin)));
    const __m128 t1 = _mm_min_ps(_mm_min_ps(_mm_max_ps(tx0, tx1), _mm_max_ps(ty0, ty1)),
                                 _mm_min_ps(_mm_max_ps(tz0, tz1), _mm_set1_ps(tmax)));
    _mm_storeu_ps(t_near, t0);
    const auto mask = static_cast<unsigned>(_mm_movemask_ps(_mm_cmple_ps(t0, t1)));
    return mask & ~empty_slots(node);
  }
};

struct Avx2Boxes {
  BSAR_TARGET_AVX2 unsigned operator()(const WideBvhNode<8>& node, const Vec3f& org,
                                       const Vec3f& inv, float tmin, float tmax,
                                       float* t_near) const {
    const __m256 ox = _mm256_set1_ps(org.x);
    const __m256 oy = _mm256_set1_ps(org.y);
    const __m256 oz = _mm256_set1_ps(org.z);
    const __m256 ix = _mm256_set1_ps(inv.x);
    const __m256 iy = _mm256_set1_ps(inv.y);
    const __m256 iz = _mm256_set1_ps(inv.z);
    const __m256 tx0 = _mm256_mul_ps(_mm256_sub_ps(_mm256_load_ps(node.lo_x), ox), ix);
    const __m256 tx1 = _mm256_mul_ps(_mm256_sub_ps(_mm256_load_ps(node.hi_x), ox), ix);
    const __m256 ty0 = _mm256_mul_ps(_mm256_sub_ps(_mm256_load_ps(node.lo_y), oy), iy);
    const __m256 ty1 = _mm256_mul_ps(_mm256_sub_ps(_mm256_load_ps(node.hi_y), oy), iy);
    const __m256 tz0 = _mm256_mul_ps(_mm256_sub_ps(_mm256_load_ps(node.lo_z), oz), iz);
    const __m256 tz1 = _mm256_mul_ps(_mm256_sub_ps(_mm256_load_ps(node.hi_z), oz), iz);
    const __m256 t0 = _mm256_max_ps(_mm256_max_ps(_mm256_min_ps(tx0, tx1), _mm256_min_ps(ty0, ty1)),
                                    _mm256_max_ps(_mm256_min_ps(tz0, tz1), _mm256_set1_ps(tmin)));
    const __m256 t1 = _mm256_min_ps(_mm256_min_ps(_mm256_max_ps(tx0, tx1), _mm256_max_ps(ty0, ty1)),
                                    _mm256_min_ps(_mm256_max_ps(tz0, tz1), _mm256_set1_ps(tmax)));
    _mm256_storeu_ps(t_near, t0);
    const auto mask = static_cast<unsigned>(_mm256_movemask_ps(_mm256_cmp_ps(t0, t1, _CMP_LE_OQ)));
    return mask & ~empty_slots(node);
  }
};
#endif

// ---------------------------------------------------------------- traversal

struct Entry {
  std::uint32_t child;
  std::uint32_t count;
  float t;
};

template <int W, typename Boxes>
inline bool traverse_closest(std::span<const WideBvhNode<W>> nodes,
                             std::span<const TriangleVerts> tris,
                             std::span<const std::uint32_t> prims, const Ray& ray, Hit& hit,
                             const Boxes& boxes) {
  const Vec3f inv = safe_inverse_direction(ray.dir);
  const WatertightRay wray(ray);
  float tmax = ray.tmax;
  bool found = false;
  std::array<Entry, kStackSize> stack;  // deliberately uninitialised
  std::size_t sp = 0;
  stack[sp++] = {0, 0, ray.tmin};

  while (sp > 0) {
    const Entry e = stack[--sp];
    if (e.t > tmax) {
      continue;
    }
    if (e.count > 0) {
      for (std::uint32_t i = e.child; i < e.child + e.count; ++i) {
        const TriangleVerts& tri = tris[i];
        float t = 0.0f;
        float u = 0.0f;
        float v = 0.0f;
        if (intersect_triangle(wray, tri.v0, tri.v1, tri.v2, ray.tmin, tmax, t, u, v)) {
          tmax = t;
          hit = {t, u, v, prims[i]};
          found = true;
        }
      }
      continue;
    }

    const WideBvhNode<W>& node = nodes[e.child];
    alignas(32) float t_near[W];
    const unsigned mask = boxes(node, ray.origin, inv, ray.tmin, tmax, t_near);
    // Push hit children far-to-near so the nearest is popped first.
    std::array<Entry, W> order;
    int k = 0;
    for (int s = 0; s < W; ++s) {
      if ((mask & (1u << s)) == 0) {
        continue;
      }
      const auto i = static_cast<std::size_t>(s);
      const Entry entry{node.child[i], node.count[i], t_near[i]};
      int j = k++;
      while (j > 0 && order[static_cast<std::size_t>(j - 1)].t < entry.t) {
        order[static_cast<std::size_t>(j)] = order[static_cast<std::size_t>(j - 1)];
        --j;
      }
      order[static_cast<std::size_t>(j)] = entry;
    }
    for (int i = 0; i < k; ++i) {
      stack[sp++] = order[static_cast<std::size_t>(i)];
    }
  }
  return found;
}

template <int W, typename Boxes>
inline bool traverse_any(std::span<const WideBvhNode<W>> nodes, std::span<const TriangleVerts> tris,
                         const Ray& ray, const Boxes& boxes) {
  const Vec3f inv = safe_inverse_direction(ray.dir);
  const WatertightRay wray(ray);
  std::array<std::uint32_t, kStackSize> stack;  // deliberately uninitialised
  std::size_t sp = 0;
  stack[sp++] = 0;
  while (sp > 0) {
    const WideBvhNode<W>& node = nodes[stack[--sp]];
    alignas(32) float t_near[W];
    const unsigned mask = boxes(node, ray.origin, inv, ray.tmin, ray.tmax, t_near);
    for (int s = 0; s < W; ++s) {
      if ((mask & (1u << s)) == 0) {
        continue;
      }
      const auto i = static_cast<std::size_t>(s);
      if (node.count[i] == 0) {
        stack[sp++] = node.child[i];
        continue;
      }
      for (std::uint32_t p = node.child[i]; p < node.child[i] + node.count[i]; ++p) {
        const TriangleVerts& tri = tris[p];
        float t = 0.0f;
        float u = 0.0f;
        float v = 0.0f;
        if (intersect_triangle(wray, tri.v0, tri.v1, tri.v2, ray.tmin, ray.tmax, t, u, v)) {
          return true;
        }
      }
    }
  }
  return false;
}

#if BSAR_X86_64
// AVX2 entry points: the whole traversal is compiled for AVX2 so the box test
// inlines into it.
BSAR_TARGET_AVX2 bool closest8_avx2(std::span<const WideBvhNode<8>> nodes,
                                    std::span<const TriangleVerts> tris,
                                    std::span<const std::uint32_t> prims, const Ray& ray,
                                    Hit& hit) {
  return traverse_closest<8>(nodes, tris, prims, ray, hit, Avx2Boxes{});
}

BSAR_TARGET_AVX2 bool any8_avx2(std::span<const WideBvhNode<8>> nodes,
                                std::span<const TriangleVerts> tris, const Ray& ray) {
  return traverse_any<8>(nodes, tris, ray, Avx2Boxes{});
}
#endif

}  // namespace

template <int W>
std::uint32_t WideBvh<W>::collapse(const Bvh& binary, std::uint32_t binary_node) {
  const auto nodes = binary.nodes();
  const auto index = static_cast<std::uint32_t>(nodes_.size());
  nodes_.push_back(empty_node<W>());

  // Gather up to W children by repeatedly opening the inner child with the
  // largest surface area.
  std::array<std::uint32_t, W> kids{};
  int n = 0;
  if (nodes[binary_node].is_leaf()) {
    kids[0] = binary_node;
    n = 1;
  } else {
    kids[0] = binary_node + 1;
    kids[1] = nodes[binary_node].index;
    n = 2;
    while (n < W) {
      int best = -1;
      float best_area = -1.0f;
      for (int i = 0; i < n; ++i) {
        const BvhNode& c = nodes[kids[static_cast<std::size_t>(i)]];
        if (!c.is_leaf() && c.bounds.surface_area() > best_area) {
          best = i;
          best_area = c.bounds.surface_area();
        }
      }
      if (best < 0) {
        break;
      }
      const std::uint32_t open = kids[static_cast<std::size_t>(best)];
      kids[static_cast<std::size_t>(best)] = open + 1;
      kids[static_cast<std::size_t>(n++)] = nodes[open].index;
    }
  }

  for (int i = 0; i < n; ++i) {
    const BvhNode& c = nodes[kids[static_cast<std::size_t>(i)]];
    if (c.is_leaf()) {
      set_slot(nodes_[index], i, c.bounds, c.index, c.count);
    } else {
      const std::uint32_t child = collapse(binary, kids[static_cast<std::size_t>(i)]);
      set_slot(nodes_[index], i, c.bounds, child, 0);  // nodes_ may have reallocated
    }
  }
  return index;
}

template <int W>
void WideBvh<W>::build(const Bvh& binary) {
  nodes_.clear();
  prims_.assign(binary.prim_indices().begin(), binary.prim_indices().end());
  tris_.assign(binary.triangles().begin(), binary.triangles().end());
  if (binary.nodes().empty()) {
    return;
  }
  nodes_.reserve(binary.nodes().size() / (W / 2) + 1);
  collapse(binary, 0);
}

template <int W>
bool WideBvh<W>::intersect(const Ray& ray, Hit& hit) const {
  if (nodes_.empty()) {
    return false;
  }
#if BSAR_X86_64
  if constexpr (W == 4) {
    return traverse_closest<4>(std::span<const Node>(nodes_), tris_, prims_, ray, hit, SseBoxes{});
  } else {
    if (cpu_has_avx2_fma()) {
      return closest8_avx2(nodes_, tris_, prims_, ray, hit);
    }
  }
#endif
  return traverse_closest<W>(std::span<const Node>(nodes_), tris_, prims_, ray, hit, ScalarBoxes{});
}

template <int W>
bool WideBvh<W>::occluded(const Ray& ray) const {
  if (nodes_.empty()) {
    return false;
  }
#if BSAR_X86_64
  if constexpr (W == 4) {
    return traverse_any<4>(std::span<const Node>(nodes_), tris_, ray, SseBoxes{});
  } else {
    if (cpu_has_avx2_fma()) {
      return any8_avx2(nodes_, tris_, ray);
    }
  }
#endif
  return traverse_any<W>(std::span<const Node>(nodes_), tris_, ray, ScalarBoxes{});
}

template class WideBvh<4>;
template class WideBvh<8>;

}  // namespace bsar
