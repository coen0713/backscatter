#include "backscatter/geometry/mesh.hpp"

#include <stdexcept>

namespace bsar {

void TriangleMesh::reserve(std::size_t vertices, std::size_t triangles) {
  vx_.reserve(vertices);
  vy_.reserve(vertices);
  vz_.reserve(vertices);
  i0_.reserve(triangles);
  i1_.reserve(triangles);
  i2_.reserve(triangles);
  material_.reserve(triangles);
}

std::uint32_t TriangleMesh::add_vertex(const Vec3f& p) {
  vx_.push_back(p.x);
  vy_.push_back(p.y);
  vz_.push_back(p.z);
  return static_cast<std::uint32_t>(vx_.size() - 1);
}

void TriangleMesh::add_triangle(std::uint32_t a, std::uint32_t b, std::uint32_t c,
                                std::uint16_t material) {
  const auto n = num_vertices();
  if (a >= n || b >= n || c >= n) {
    throw std::out_of_range("TriangleMesh::add_triangle: vertex index out of range");
  }
  i0_.push_back(a);
  i1_.push_back(b);
  i2_.push_back(c);
  material_.push_back(material);
}

void TriangleMesh::append(const TriangleMesh& other) {
  const auto offset = static_cast<std::uint32_t>(num_vertices());
  vx_.insert(vx_.end(), other.vx_.begin(), other.vx_.end());
  vy_.insert(vy_.end(), other.vy_.begin(), other.vy_.end());
  vz_.insert(vz_.end(), other.vz_.begin(), other.vz_.end());
  for (std::size_t t = 0; t < other.num_triangles(); ++t) {
    i0_.push_back(other.i0_[t] + offset);
    i1_.push_back(other.i1_[t] + offset);
    i2_.push_back(other.i2_[t] + offset);
  }
  material_.insert(material_.end(), other.material_.begin(), other.material_.end());
}

Vec3f TriangleMesh::normal(std::size_t t) const {
  const auto [a, b, c] = triangle(t);
  return normalize(cross(b - a, c - a));
}

float TriangleMesh::area(std::size_t t) const {
  const auto [a, b, c] = triangle(t);
  return 0.5f * length(cross(b - a, c - a));
}

Aabb TriangleMesh::triangle_bounds(std::size_t t) const {
  const auto [a, b, c] = triangle(t);
  Aabb box;
  box.expand(a);
  box.expand(b);
  box.expand(c);
  return box;
}

Aabb TriangleMesh::bounds() const {
  Aabb box;
  for (std::size_t i = 0; i < num_vertices(); ++i) {
    box.expand(Vec3f{vx_[i], vy_[i], vz_[i]});
  }
  return box;
}

}  // namespace bsar
