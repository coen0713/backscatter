#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "backscatter/math/aabb.hpp"
#include "backscatter/math/vec3.hpp"

namespace bsar {

/// Indexed triangle mesh in structure-of-arrays layout. Coordinates are
/// metres in the scene-local east-north-up frame.
class TriangleMesh {
 public:
  [[nodiscard]] std::size_t num_vertices() const { return vx_.size(); }
  [[nodiscard]] std::size_t num_triangles() const { return i0_.size(); }
  [[nodiscard]] bool empty() const { return i0_.empty(); }

  void reserve(std::size_t vertices, std::size_t triangles);

  std::uint32_t add_vertex(const Vec3f& p);
  void add_triangle(std::uint32_t a, std::uint32_t b, std::uint32_t c, std::uint16_t material);
  void append(const TriangleMesh& other);

  [[nodiscard]] Vec3f vertex(std::uint32_t i) const { return {vx_[i], vy_[i], vz_[i]}; }
  [[nodiscard]] std::array<Vec3f, 3> triangle(std::size_t t) const {
    return {vertex(i0_[t]), vertex(i1_[t]), vertex(i2_[t])};
  }
  [[nodiscard]] std::uint16_t material(std::size_t t) const { return material_[t]; }

  /// Unit geometric normal following the counter-clockwise winding rule.
  [[nodiscard]] Vec3f normal(std::size_t t) const;
  [[nodiscard]] float area(std::size_t t) const;
  [[nodiscard]] Aabb triangle_bounds(std::size_t t) const;
  [[nodiscard]] Aabb bounds() const;

 private:
  std::vector<float> vx_, vy_, vz_;
  std::vector<std::uint32_t> i0_, i1_, i2_;
  std::vector<std::uint16_t> material_;
};

}  // namespace bsar
