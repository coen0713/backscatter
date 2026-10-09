#include "backscatter/scene/buildings.hpp"

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace bsar {
namespace {

double cross2(const Vec2d& o, const Vec2d& a, const Vec2d& b) {
  return (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x);
}

bool point_in_triangle(const Vec2d& p, const Vec2d& a, const Vec2d& b, const Vec2d& c) {
  return cross2(a, b, p) >= 0.0 && cross2(b, c, p) >= 0.0 && cross2(c, a, p) >= 0.0;
}

std::vector<Vec2d> clean_ring(std::span<const Vec2d> ring) {
  std::vector<Vec2d> out;
  for (const Vec2d& p : ring) {
    if (out.empty() || p.x != out.back().x || p.y != out.back().y) {
      out.push_back(p);
    }
  }
  while (out.size() > 1 && out.front().x == out.back().x && out.front().y == out.back().y) {
    out.pop_back();
  }
  return out;
}

}  // namespace

double signed_area(std::span<const Vec2d> ring) {
  double a = 0.0;
  for (std::size_t i = 0; i < ring.size(); ++i) {
    const Vec2d& p = ring[i];
    const Vec2d& q = ring[(i + 1) % ring.size()];
    a += p.x * q.y - q.x * p.y;
  }
  return 0.5 * a;
}

std::vector<std::uint32_t> triangulate_polygon(std::span<const Vec2d> ring) {
  std::vector<std::uint32_t> out;
  const std::size_t n = ring.size();
  if (n < 3) {
    return out;
  }
  std::vector<std::uint32_t> idx(n);
  for (std::size_t i = 0; i < n; ++i) {
    idx[i] = static_cast<std::uint32_t>(i);
  }
  if (signed_area(ring) < 0.0) {
    std::reverse(idx.begin(), idx.end());
  }

  std::size_t guard = 0;
  while (idx.size() > 3 && guard < 4 * n * n) {
    ++guard;
    bool clipped = false;
    for (std::size_t k = 0; k < idx.size(); ++k) {
      const std::uint32_t ia = idx[(k + idx.size() - 1) % idx.size()];
      const std::uint32_t ib = idx[k];
      const std::uint32_t ic = idx[(k + 1) % idx.size()];
      const Vec2d& a = ring[ia];
      const Vec2d& b = ring[ib];
      const Vec2d& c = ring[ic];
      if (cross2(a, b, c) <= 0.0) {
        continue;  // reflex or degenerate corner
      }
      bool ear = true;
      for (const std::uint32_t j : idx) {
        if (j != ia && j != ib && j != ic && point_in_triangle(ring[j], a, b, c)) {
          ear = false;
          break;
        }
      }
      if (ear) {
        out.insert(out.end(), {ia, ib, ic});
        idx.erase(idx.begin() + static_cast<std::ptrdiff_t>(k));
        clipped = true;
        break;
      }
    }
    if (!clipped) {
      break;  // self-intersecting input; fall back to a fan for the rest
    }
  }
  for (std::size_t k = 1; k + 1 < idx.size(); ++k) {
    out.insert(out.end(), {idx[0], idx[k], idx[k + 1]});
  }
  return out;
}

void extrude_footprint(std::span<const Vec2d> ring_in, double base, double top,
                       std::uint16_t material, TriangleMesh& out) {
  std::vector<Vec2d> ring = clean_ring(ring_in);
  if (ring.size() < 3 || top <= base) {
    return;
  }
  if (signed_area(ring) < 0.0) {
    std::reverse(ring.begin(), ring.end());
  }
  const std::size_t n = ring.size();
  const auto first = static_cast<std::uint32_t>(out.num_vertices());
  for (const Vec2d& p : ring) {
    out.add_vertex(Vec3f(Vec3d{p.x, p.y, base}));
  }
  for (const Vec2d& p : ring) {
    out.add_vertex(Vec3f(Vec3d{p.x, p.y, top}));
  }
  // Walls: counter-clockwise ring => these windings face outwards.
  for (std::size_t i = 0; i < n; ++i) {
    const auto a = first + static_cast<std::uint32_t>(i);
    const auto b = first + static_cast<std::uint32_t>((i + 1) % n);
    const auto at = a + static_cast<std::uint32_t>(n);
    const auto bt = b + static_cast<std::uint32_t>(n);
    out.add_triangle(a, b, bt, material);
    out.add_triangle(a, bt, at, material);
  }
  // Roof.
  const auto tris = triangulate_polygon(ring);
  const auto roof = first + static_cast<std::uint32_t>(n);
  for (std::size_t k = 0; k + 2 < tris.size(); k += 3) {
    out.add_triangle(roof + tris[k], roof + tris[k + 1], roof + tris[k + 2], material);
  }
}

FootprintSet load_footprints(const std::filesystem::path& path) {
  std::ifstream in(path);
  if (!in) {
    throw std::runtime_error("load_footprints: cannot open " + path.string());
  }
  FootprintSet set;
  std::string line;
  Footprint current;
  bool open = false;
  std::size_t line_no = 0;
  while (std::getline(in, line)) {
    ++line_no;
    if (const auto hash = line.find('#'); hash != std::string::npos) {
      line.erase(hash);
    }
    std::istringstream ls(line);
    std::string word;
    if (!(ls >> word)) {
      continue;
    }
    if (word == "crs") {
      std::string crs;
      ls >> crs;
      if (crs != "wgs84" && crs != "enu") {
        throw std::runtime_error("load_footprints: unknown crs '" + crs + "'");
      }
      set.geographic = crs == "wgs84";
    } else if (word == "building") {
      current = Footprint{};
      if (!(ls >> current.height)) {
        throw std::runtime_error("load_footprints: missing height on line " +
                                 std::to_string(line_no));
      }
      std::string material;
      if (ls >> material) {
        current.material = material;
      }
      open = true;
    } else if (word == "end") {
      if (open) {
        set.footprints.push_back(std::exchange(current, Footprint{}));
      }
      open = false;
    } else if (open) {
      Vec2d p;
      p.x = std::stod(word);
      if (!(ls >> p.y)) {
        throw std::runtime_error("load_footprints: bad vertex on line " + std::to_string(line_no));
      }
      current.ring.push_back(p);
    } else {
      throw std::runtime_error("load_footprints: unexpected '" + word + "' on line " +
                               std::to_string(line_no));
    }
  }
  return set;
}

}  // namespace bsar
