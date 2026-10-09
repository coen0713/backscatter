#pragma once

#include <cstdint>
#include <span>
#include <vector>

#include "backscatter/geometry/bvh.hpp"
#include "backscatter/geometry/ray.hpp"

namespace bsar {

/// W-wide node: the boxes of all W children are stored lane-wise so one SIMD
/// slab test checks them together.
template <int W>
struct alignas(64) WideBvhNode {
  static constexpr std::uint32_t kEmpty = 0xFFFFFFFFu;

  float lo_x[W], hi_x[W];
  float lo_y[W], hi_y[W];
  float lo_z[W], hi_z[W];
  std::uint32_t child[W];  // inner: node index; leaf: first triangle; kEmpty: unused slot
  std::uint32_t count[W];  // 0 for inner children, triangle count for leaves
};

/// W-wide BVH obtained by collapsing a binary SAH BVH: each node absorbs
/// grandchildren (largest surface area first) until it has W children.
///
/// * W = 4: one SSE slab test per node (SSE2 is baseline on x86-64).
/// * W = 8: one AVX2 slab test per node, used when the CPU supports AVX2
///   (checked at run time); otherwise the same traversal runs with a scalar
///   box loop.
///
/// Every variant reports the same closest hit as the binary BVH.
template <int W>
class WideBvh {
 public:
  static_assert(W == 4 || W == 8, "WideBvh supports widths 4 and 8");
  using Node = WideBvhNode<W>;

  void build(const Bvh& binary);

  bool intersect(const Ray& ray, Hit& hit) const;
  [[nodiscard]] bool occluded(const Ray& ray) const;

  [[nodiscard]] std::span<const Node> nodes() const { return nodes_; }

 private:
  std::uint32_t collapse(const Bvh& binary, std::uint32_t binary_node);

  std::vector<Node> nodes_;
  std::vector<std::uint32_t> prims_;
  std::vector<TriangleVerts> tris_;
};

using Bvh4 = WideBvh<4>;
using Bvh8 = WideBvh<8>;

extern template class WideBvh<4>;
extern template class WideBvh<8>;

}  // namespace bsar
