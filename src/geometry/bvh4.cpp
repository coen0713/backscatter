#include "backscatter/geometry/bvh4.hpp"

#include <algorithm>
#include <array>
#include <limits>

#include "backscatter/geometry/ray_box.hpp"
#include "backscatter/geometry/triangle.hpp"

#if defined(__SSE2__) || defined(_M_X64) || (defined(_M_IX86_FP) && _M_IX86_FP >= 2)
#define BSAR_HAVE_SSE 1
#include <emmintrin.h>
#else
#define BSAR_HAVE_SSE 0
#endif

namespace bsar {
namespace {

constexpr std::size_t kStackSize = 256;

void set_slot(Bvh4Node& node, int slot, const Aabb& box, std::uint32_t child, std::uint32_t count) {
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

Bvh4Node empty_node() {
  Bvh4Node node{};
  for (int s = 0; s < 4; ++s) {
    // Inverted boxes never pass the slab test; the slot mask below also
    // excludes them explicitly.
    set_slot(node, s, Aabb{}, Bvh4Node::kEmpty, 0);
  }
  return node;
}

/// Slab test of one ray against the four child boxes. Returns a bit mask of
/// the children hit, and their entry distances.
inline unsigned intersect_children(const Bvh4Node& node, const Vec3f& org, const Vec3f& inv,
                                   float tmin, float tmax, std::array<float, 4>& t_near) {
#if BSAR_HAVE_SSE
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
  _mm_storeu_ps(t_near.data(), t0);
  auto mask = static_cast<unsigned>(_mm_movemask_ps(_mm_cmple_ps(t0, t1)));
#else
  unsigned mask = 0;
  for (int s = 0; s < 4; ++s) {
    const auto i = static_cast<std::size_t>(s);
    Aabb box;
    box.lo = {node.lo_x[i], node.lo_y[i], node.lo_z[i]};
    box.hi = {node.hi_x[i], node.hi_y[i], node.hi_z[i]};
    if (intersect_aabb(box, org, inv, tmin, tmax, t_near[i])) {
      mask |= 1u << s;
    }
  }
#endif
  for (int s = 0; s < 4; ++s) {
    if (node.child[static_cast<std::size_t>(s)] == Bvh4Node::kEmpty) {
      mask &= ~(1u << s);
    }
  }
  return mask;
}

}  // namespace

std::uint32_t Bvh4::collapse(const Bvh& binary, std::uint32_t binary_node) {
  const auto nodes = binary.nodes();
  const auto index = static_cast<std::uint32_t>(nodes_.size());
  nodes_.push_back(empty_node());

  // Gather up to four children by repeatedly opening the inner child with the
  // largest surface area.
  std::array<std::uint32_t, 4> kids{};
  int n = 0;
  if (nodes[binary_node].is_leaf()) {
    kids[0] = binary_node;
    n = 1;
  } else {
    kids[0] = binary_node + 1;
    kids[1] = nodes[binary_node].index;
    n = 2;
    while (n < 4) {
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

void Bvh4::build(const Bvh& binary) {
  nodes_.clear();
  prims_.assign(binary.prim_indices().begin(), binary.prim_indices().end());
  tris_.assign(binary.triangles().begin(), binary.triangles().end());
  if (binary.nodes().empty()) {
    return;
  }
  nodes_.reserve(binary.nodes().size() / 2 + 1);
  collapse(binary, 0);
}

bool Bvh4::intersect(const Ray& ray, Hit& hit) const {
  if (nodes_.empty()) {
    return false;
  }
  const Vec3f inv = safe_inverse_direction(ray.dir);
  const WatertightRay wray(ray);
  float tmax = ray.tmax;
  bool found = false;

  struct Entry {
    std::uint32_t child;
    std::uint32_t count;
    float t;
  };
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
        const TriangleVerts& tri = tris_[i];
        float t = 0.0f;
        float u = 0.0f;
        float v = 0.0f;
        if (intersect_triangle(wray, tri.v0, tri.v1, tri.v2, ray.tmin, tmax, t, u, v)) {
          tmax = t;
          hit = {t, u, v, prims_[i]};
          found = true;
        }
      }
      continue;
    }

    const Bvh4Node& node = nodes_[e.child];
    std::array<float, 4> t_near{};
    const unsigned mask = intersect_children(node, ray.origin, inv, ray.tmin, tmax, t_near);
    // Push hit children far-to-near so the nearest is popped first.
    std::array<Entry, 4> hits{};
    int k = 0;
    for (int s = 0; s < 4; ++s) {
      if ((mask & (1u << s)) != 0) {
        const auto i = static_cast<std::size_t>(s);
        Entry entry{node.child[i], node.count[i], t_near[i]};
        int j = k++;
        while (j > 0 && hits[static_cast<std::size_t>(j - 1)].t < entry.t) {
          hits[static_cast<std::size_t>(j)] = hits[static_cast<std::size_t>(j - 1)];
          --j;
        }
        hits[static_cast<std::size_t>(j)] = entry;
      }
    }
    for (int i = 0; i < k; ++i) {
      stack[sp++] = hits[static_cast<std::size_t>(i)];
    }
  }
  return found;
}

bool Bvh4::occluded(const Ray& ray) const {
  if (nodes_.empty()) {
    return false;
  }
  const Vec3f inv = safe_inverse_direction(ray.dir);
  const WatertightRay wray(ray);
  std::array<std::uint32_t, kStackSize> stack;  // deliberately uninitialised
  std::size_t sp = 0;
  stack[sp++] = 0;
  while (sp > 0) {
    const Bvh4Node& node = nodes_[stack[--sp]];
    std::array<float, 4> t_near{};
    const unsigned mask = intersect_children(node, ray.origin, inv, ray.tmin, ray.tmax, t_near);
    for (int s = 0; s < 4; ++s) {
      if ((mask & (1u << s)) == 0) {
        continue;
      }
      const auto i = static_cast<std::size_t>(s);
      if (node.count[i] == 0) {
        stack[sp++] = node.child[i];
        continue;
      }
      for (std::uint32_t p = node.child[i]; p < node.child[i] + node.count[i]; ++p) {
        const TriangleVerts& tri = tris_[p];
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

}  // namespace bsar
