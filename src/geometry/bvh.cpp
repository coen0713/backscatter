#include "backscatter/geometry/bvh.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <limits>
#include <utility>

#include "backscatter/geometry/ray_box.hpp"
#include "backscatter/geometry/triangle.hpp"

namespace bsar {
namespace {

constexpr std::size_t kStackSize = 192;
// Beyond this depth the builder switches to balanced splits, which bounds the
// tree depth (and so the traversal stack) for any input.
constexpr std::size_t kBalancedDepth = 96;
constexpr int kMaxBins = 64;

struct BuildPrim {
  Aabb box;
  Vec3f centroid;
};

class Builder {
 public:
  Builder(const BvhBuildConfig& config, std::vector<BuildPrim> prims,
          std::vector<std::uint32_t>& order, std::vector<BvhNode>& nodes)
      : config_(config), prims_(std::move(prims)), order_(order), nodes_(nodes) {}

  void run() {
    nodes_.reserve(
        2 * prims_.size() / static_cast<std::size_t>(std::max(1, config_.max_leaf_size)) + 1);
    build(0, static_cast<std::uint32_t>(order_.size()), 0);
  }

  [[nodiscard]] std::size_t max_depth() const { return max_depth_; }

 private:
  struct Split {
    int axis = -1;
    int bin = 0;
    float cost = std::numeric_limits<float>::infinity();
  };

  std::uint32_t make_leaf(std::uint32_t node, std::uint32_t begin, std::uint32_t end) {
    nodes_[node].index = begin;
    nodes_[node].count = static_cast<std::uint16_t>(end - begin);
    return node;
  }

  Split find_sah_split(std::uint32_t begin, std::uint32_t end, const Aabb& bounds,
                       const Aabb& centroid_bounds) const {
    const int num_bins = std::clamp(config_.num_bins, 2, kMaxBins);
    const float parent_area = bounds.surface_area();
    Split best;
    for (int axis = 0; axis < 3; ++axis) {
      const auto ax = static_cast<std::size_t>(axis);
      const float lo = centroid_bounds.lo[ax];
      const float extent = centroid_bounds.hi[ax] - lo;
      if (!(extent > 0.0f)) {
        continue;
      }
      std::array<Aabb, kMaxBins> bin_box{};
      std::array<std::uint32_t, kMaxBins> bin_count{};
      const float scale = static_cast<float>(num_bins) / extent;
      for (std::uint32_t i = begin; i < end; ++i) {
        const BuildPrim& p = prims_[order_[i]];
        const int b = std::min(num_bins - 1, static_cast<int>((p.centroid[ax] - lo) * scale));
        bin_box[static_cast<std::size_t>(b)].expand(p.box);
        ++bin_count[static_cast<std::size_t>(b)];
      }
      // Sweep from the right to get the area/count of every right partition.
      std::array<float, kMaxBins> right_area{};
      std::array<std::uint32_t, kMaxBins> right_count{};
      Aabb acc;
      std::uint32_t n = 0;
      for (int b = num_bins - 1; b > 0; --b) {
        acc.expand(bin_box[static_cast<std::size_t>(b)]);
        n += bin_count[static_cast<std::size_t>(b)];
        right_area[static_cast<std::size_t>(b - 1)] = acc.surface_area();
        right_count[static_cast<std::size_t>(b - 1)] = n;
      }
      acc = Aabb{};
      n = 0;
      for (int b = 0; b < num_bins - 1; ++b) {
        acc.expand(bin_box[static_cast<std::size_t>(b)]);
        n += bin_count[static_cast<std::size_t>(b)];
        const std::uint32_t nr = right_count[static_cast<std::size_t>(b)];
        if (n == 0 || nr == 0) {
          continue;
        }
        const float cost = config_.traversal_cost +
                           config_.intersection_cost *
                               (acc.surface_area() * static_cast<float>(n) +
                                right_area[static_cast<std::size_t>(b)] * static_cast<float>(nr)) /
                               parent_area;
        if (cost < best.cost) {
          best = {axis, b, cost};
        }
      }
    }
    return best;
  }

  std::uint32_t split_balanced(std::uint32_t begin, std::uint32_t end, int axis) {
    const std::uint32_t mid = begin + (end - begin) / 2;
    const auto ax = static_cast<std::size_t>(axis);
    std::nth_element(order_.begin() + begin, order_.begin() + mid, order_.begin() + end,
                     [&](std::uint32_t a, std::uint32_t b) {
                       return prims_[a].centroid[ax] < prims_[b].centroid[ax];
                     });
    return mid;
  }

  std::uint32_t build(std::uint32_t begin, std::uint32_t end, std::size_t depth) {
    max_depth_ = std::max(max_depth_, depth);
    const auto node = static_cast<std::uint32_t>(nodes_.size());
    nodes_.emplace_back();

    Aabb bounds;
    Aabb centroid_bounds;
    for (std::uint32_t i = begin; i < end; ++i) {
      bounds.expand(prims_[order_[i]].box);
      centroid_bounds.expand(prims_[order_[i]].centroid);
    }
    nodes_[node].bounds = bounds;

    const std::uint32_t count = end - begin;
    const auto max_leaf = static_cast<std::uint32_t>(std::max(1, config_.max_leaf_size));
    if (count <= 1) {
      return make_leaf(node, begin, end);
    }

    int axis = centroid_bounds.largest_axis();
    std::uint32_t mid = begin;
    const bool degenerate = !(centroid_bounds.extent()[static_cast<std::size_t>(axis)] > 0.0f);

    if (degenerate || depth >= kBalancedDepth ||
        config_.strategy == BvhBuildStrategy::ObjectMedian) {
      if (count <= max_leaf) {
        return make_leaf(node, begin, end);
      }
      mid = split_balanced(begin, end, axis);
    } else {
      const Split split = find_sah_split(begin, end, bounds, centroid_bounds);
      const float leaf_cost = config_.intersection_cost * static_cast<float>(count);
      if (count <= max_leaf && (split.axis < 0 || split.cost >= leaf_cost)) {
        return make_leaf(node, begin, end);
      }
      if (split.axis < 0) {
        mid = split_balanced(begin, end, axis);
      } else {
        axis = split.axis;
        const auto ax = static_cast<std::size_t>(axis);
        const int num_bins = std::clamp(config_.num_bins, 2, kMaxBins);
        const float lo = centroid_bounds.lo[ax];
        const float scale = static_cast<float>(num_bins) / (centroid_bounds.hi[ax] - lo);
        const auto it =
            std::partition(order_.begin() + begin, order_.begin() + end, [&](std::uint32_t p) {
              const int b =
                  std::min(num_bins - 1, static_cast<int>((prims_[p].centroid[ax] - lo) * scale));
              return b <= split.bin;
            });
        mid = static_cast<std::uint32_t>(it - order_.begin());
        if (mid == begin || mid == end) {
          mid = split_balanced(begin, end, axis);
        }
      }
    }

    build(begin, mid, depth + 1);  // left child lands at node + 1
    const std::uint32_t right = build(mid, end, depth + 1);
    nodes_[node].index = right;
    nodes_[node].axis = static_cast<std::uint16_t>(axis);
    return node;
  }

  const BvhBuildConfig& config_;
  std::vector<BuildPrim> prims_;
  std::vector<std::uint32_t>& order_;
  std::vector<BvhNode>& nodes_;
  std::size_t max_depth_ = 0;
};

}  // namespace

double compute_sah_cost(std::span<const BvhNode> nodes, float traversal_cost,
                        float intersection_cost) {
  if (nodes.empty()) {
    return 0.0;
  }
  const double root_area = nodes.front().bounds.surface_area();
  if (root_area <= 0.0) {
    return 0.0;
  }
  double cost = 0.0;
  for (const BvhNode& n : nodes) {
    const double rel = n.bounds.surface_area() / root_area;
    cost += n.is_leaf() ? rel * intersection_cost * n.count : rel * traversal_cost;
  }
  return cost;
}

void Bvh::build(const TriangleMesh& mesh, const BvhBuildConfig& config) {
  const auto start = std::chrono::steady_clock::now();
  nodes_.clear();
  prims_.clear();
  tris_.clear();
  stats_ = {};

  const std::size_t n = mesh.num_triangles();
  if (n == 0) {
    return;
  }

  std::vector<BuildPrim> prims(n);
  for (std::size_t t = 0; t < n; ++t) {
    prims[t].box = mesh.triangle_bounds(t);
    prims[t].centroid = prims[t].box.center();
  }
  prims_.resize(n);
  for (std::size_t t = 0; t < n; ++t) {
    prims_[t] = static_cast<std::uint32_t>(t);
  }

  Builder builder(config, std::move(prims), prims_, nodes_);
  builder.run();

  tris_.resize(n);
  for (std::size_t i = 0; i < n; ++i) {
    const auto [a, b, c] = mesh.triangle(prims_[i]);
    tris_[i] = {a, b, c};
  }

  stats_.num_nodes = nodes_.size();
  stats_.num_leaves = static_cast<std::size_t>(
      std::count_if(nodes_.begin(), nodes_.end(), [](const BvhNode& nd) { return nd.is_leaf(); }));
  stats_.max_depth = builder.max_depth();
  stats_.sah_cost = compute_sah_cost(nodes_, config.traversal_cost, config.intersection_cost);
  stats_.build_seconds =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}

bool Bvh::intersect(const Ray& ray, Hit& hit) const {
  if (nodes_.empty()) {
    return false;
  }
  const Vec3f inv = safe_inverse_direction(ray.dir);
  const WatertightRay wray(ray);
  float tmax = ray.tmax;
  float t_entry = 0.0f;
  if (!intersect_aabb(nodes_[0].bounds, ray.origin, inv, ray.tmin, tmax, t_entry)) {
    return false;
  }

  struct StackEntry {
    std::uint32_t node;
    float t;
  };
  std::array<StackEntry, kStackSize> stack;  // deliberately uninitialised
  std::size_t sp = 0;
  std::uint32_t current = 0;
  bool found = false;

  for (;;) {
    const BvhNode& node = nodes_[current];
    if (node.is_leaf()) {
      for (std::uint32_t i = node.index; i < node.index + node.count; ++i) {
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
    } else {
      std::uint32_t near = current + 1;
      std::uint32_t far = node.index;
      float t_near = 0.0f;
      float t_far = 0.0f;
      bool hit_near = intersect_aabb(nodes_[near].bounds, ray.origin, inv, ray.tmin, tmax, t_near);
      bool hit_far = intersect_aabb(nodes_[far].bounds, ray.origin, inv, ray.tmin, tmax, t_far);
      if (hit_near && hit_far) {
        if (t_far < t_near) {
          std::swap(near, far);
          std::swap(t_near, t_far);
        }
        stack[sp++] = {far, t_far};
        current = near;
        continue;
      }
      if (hit_near || hit_far) {
        current = hit_near ? near : far;
        continue;
      }
    }
    // Pop, skipping nodes that are now farther than the closest hit.
    for (;;) {
      if (sp == 0) {
        return found;
      }
      const auto [idx, t] = stack[--sp];
      if (t <= tmax) {
        current = idx;
        break;
      }
    }
  }
}

bool Bvh::occluded(const Ray& ray) const {
  if (nodes_.empty()) {
    return false;
  }
  const Vec3f inv = safe_inverse_direction(ray.dir);
  const WatertightRay wray(ray);
  std::array<std::uint32_t, kStackSize> stack;  // deliberately uninitialised
  std::size_t sp = 0;
  stack[sp++] = 0;
  while (sp > 0) {
    const BvhNode& node = nodes_[stack[--sp]];
    float t_entry = 0.0f;
    if (!intersect_aabb(node.bounds, ray.origin, inv, ray.tmin, ray.tmax, t_entry)) {
      continue;
    }
    if (node.is_leaf()) {
      for (std::uint32_t i = node.index; i < node.index + node.count; ++i) {
        const TriangleVerts& tri = tris_[i];
        float t = 0.0f;
        float u = 0.0f;
        float v = 0.0f;
        if (intersect_triangle(wray, tri.v0, tri.v1, tri.v2, ray.tmin, ray.tmax, t, u, v)) {
          return true;
        }
      }
    } else {
      const auto self = static_cast<std::uint32_t>(&node - nodes_.data());
      stack[sp++] = node.index;
      stack[sp++] = self + 1;
    }
  }
  return false;
}

}  // namespace bsar
