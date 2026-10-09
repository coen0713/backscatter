#pragma once

#include <cstdint>
#include <span>
#include <vector>

#include "backscatter/geometry/bvh.hpp"
#include "backscatter/geometry/ray.hpp"

namespace bsar {

/// Four-wide node: the boxes of all four children are stored lane-wise so one
/// SIMD slab test checks them together.
struct alignas(64) Bvh4Node {
  static constexpr std::uint32_t kEmpty = 0xFFFFFFFFu;

  float lo_x[4], hi_x[4];
  float lo_y[4], hi_y[4];
  float lo_z[4], hi_z[4];
  std::uint32_t child[4];  // inner: node index; leaf: first triangle; kEmpty: unused slot
  std::uint32_t count[4];  // 0 for inner children, triangle count for leaves
};

/// Four-wide BVH obtained by collapsing a binary SAH BVH: each node absorbs
/// grandchildren (largest surface area first) until it has four children.
/// Traversal uses SSE when available and a scalar loop otherwise; both give
/// identical hits.
class Bvh4 {
 public:
  void build(const Bvh& binary);

  bool intersect(const Ray& ray, Hit& hit) const;
  [[nodiscard]] bool occluded(const Ray& ray) const;

  [[nodiscard]] std::span<const Bvh4Node> nodes() const { return nodes_; }

 private:
  std::uint32_t collapse(const Bvh& binary, std::uint32_t binary_node);

  std::vector<Bvh4Node> nodes_;
  std::vector<std::uint32_t> prims_;
  std::vector<TriangleVerts> tris_;
};

}  // namespace bsar
