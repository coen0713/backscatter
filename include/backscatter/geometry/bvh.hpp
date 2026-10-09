#pragma once

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "backscatter/geometry/mesh.hpp"
#include "backscatter/geometry/ray.hpp"
#include "backscatter/math/aabb.hpp"

namespace bsar {

enum class BvhBuildStrategy {
  BinnedSah,     // surface area heuristic evaluated over equal-width centroid bins
  ObjectMedian,  // split at the median centroid of the largest axis (baseline)
};

struct BvhBuildConfig {
  BvhBuildStrategy strategy = BvhBuildStrategy::BinnedSah;
  int num_bins = 16;
  int max_leaf_size = 4;  // leaves never exceed this many triangles
  float traversal_cost = 1.0f;
  float intersection_cost = 1.0f;
};

/// 32-byte node. Inner nodes store their left child immediately after
/// themselves (depth-first layout) and the right child index in `index`.
struct BvhNode {
  Aabb bounds;
  std::uint32_t index = 0;  // leaf: first triangle; inner: right child
  std::uint16_t count = 0;  // leaf: triangle count; inner: 0
  std::uint16_t axis = 0;   // inner: split axis, used for ordered traversal

  [[nodiscard]] bool is_leaf() const { return count > 0; }
};

struct BvhStats {
  std::size_t num_nodes = 0;
  std::size_t num_leaves = 0;
  std::size_t max_depth = 0;
  double sah_cost = 0.0;  // expected cost per ray, normalised by root area
  double build_seconds = 0.0;
};

/// Triangle vertices stored in leaf order for cache-friendly traversal.
struct TriangleVerts {
  Vec3f v0, v1, v2;
};

/// Binary bounding volume hierarchy over a triangle mesh.
class Bvh {
 public:
  void build(const TriangleMesh& mesh, const BvhBuildConfig& config = {});

  /// Closest hit in (ray.tmin, ray.tmax). `hit.prim` is the mesh triangle index.
  bool intersect(const Ray& ray, Hit& hit) const;
  /// True if any triangle is hit in (ray.tmin, ray.tmax).
  [[nodiscard]] bool occluded(const Ray& ray) const;

  [[nodiscard]] std::span<const BvhNode> nodes() const { return nodes_; }
  [[nodiscard]] std::span<const std::uint32_t> prim_indices() const { return prims_; }
  [[nodiscard]] std::span<const TriangleVerts> triangles() const { return tris_; }
  [[nodiscard]] const BvhStats& stats() const { return stats_; }
  [[nodiscard]] Aabb bounds() const { return nodes_.empty() ? Aabb{} : nodes_.front().bounds; }

 private:
  std::vector<BvhNode> nodes_;
  std::vector<std::uint32_t> prims_;  // leaf order -> mesh triangle index
  std::vector<TriangleVerts> tris_;   // leaf order
  BvhStats stats_;
};

/// Expected traversal cost of a built hierarchy under the SAH cost model.
double compute_sah_cost(std::span<const BvhNode> nodes, float traversal_cost,
                        float intersection_cost);

}  // namespace bsar
