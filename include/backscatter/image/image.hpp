#pragma once

#include <cstddef>
#include <vector>

namespace bsar {

/// Row-major 2D array. For SAR products x is range (columns) and y is
/// azimuth (rows).
template <typename T>
struct Image {
  std::size_t width = 0;
  std::size_t height = 0;
  std::vector<T> data;

  Image() = default;
  Image(std::size_t w, std::size_t h, T fill = T{}) : width(w), height(h), data(w * h, fill) {}

  [[nodiscard]] T& at(std::size_t x, std::size_t y) { return data[y * width + x]; }
  [[nodiscard]] const T& at(std::size_t x, std::size_t y) const { return data[y * width + x]; }
  [[nodiscard]] T* row(std::size_t y) { return data.data() + y * width; }
  [[nodiscard]] const T* row(std::size_t y) const { return data.data() + y * width; }
  [[nodiscard]] std::size_t size() const { return data.size(); }
};

}  // namespace bsar
