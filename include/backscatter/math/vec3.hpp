#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace bsar {

/// Minimal 3-vector. `float` is used for scene geometry and BVH traversal,
/// `double` for sensor geometry, where path lengths of ~10^6 m must be
/// resolved to a fraction of a wavelength.
template <typename T>
struct Vec3 {
  T x{};
  T y{};
  T z{};

  constexpr Vec3() = default;
  constexpr Vec3(T x_, T y_, T z_) : x(x_), y(y_), z(z_) {}

  template <typename U>
  constexpr explicit Vec3(const Vec3<U>& o)
      : x(static_cast<T>(o.x)), y(static_cast<T>(o.y)), z(static_cast<T>(o.z)) {}

  constexpr T operator[](std::size_t i) const { return i == 0 ? x : (i == 1 ? y : z); }
  constexpr T& operator[](std::size_t i) { return i == 0 ? x : (i == 1 ? y : z); }

  constexpr Vec3& operator+=(const Vec3& o) {
    x += o.x;
    y += o.y;
    z += o.z;
    return *this;
  }
  constexpr Vec3& operator-=(const Vec3& o) {
    x -= o.x;
    y -= o.y;
    z -= o.z;
    return *this;
  }
  constexpr Vec3& operator*=(T s) {
    x *= s;
    y *= s;
    z *= s;
    return *this;
  }
  constexpr Vec3& operator/=(T s) {
    x /= s;
    y /= s;
    z /= s;
    return *this;
  }
  constexpr bool operator==(const Vec3&) const = default;
};

using Vec3f = Vec3<float>;
using Vec3d = Vec3<double>;

template <typename T>
constexpr Vec3<T> operator+(Vec3<T> a, const Vec3<T>& b) {
  return a += b;
}
template <typename T>
constexpr Vec3<T> operator-(Vec3<T> a, const Vec3<T>& b) {
  return a -= b;
}
template <typename T>
constexpr Vec3<T> operator-(const Vec3<T>& a) {
  return {-a.x, -a.y, -a.z};
}
template <typename T>
constexpr Vec3<T> operator*(Vec3<T> a, T s) {
  return a *= s;
}
template <typename T>
constexpr Vec3<T> operator*(T s, Vec3<T> a) {
  return a *= s;
}
template <typename T>
constexpr Vec3<T> operator/(Vec3<T> a, T s) {
  return a /= s;
}

template <typename T>
constexpr T dot(const Vec3<T>& a, const Vec3<T>& b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

template <typename T>
constexpr Vec3<T> cross(const Vec3<T>& a, const Vec3<T>& b) {
  return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

template <typename T>
constexpr Vec3<T> hadamard(const Vec3<T>& a, const Vec3<T>& b) {
  return {a.x * b.x, a.y * b.y, a.z * b.z};
}

template <typename T>
constexpr T length_squared(const Vec3<T>& a) {
  return dot(a, a);
}

template <typename T>
T length(const Vec3<T>& a) {
  return std::sqrt(dot(a, a));
}

template <typename T>
Vec3<T> normalize(const Vec3<T>& a) {
  const T len = length(a);
  return len > T(0) ? a / len : a;
}

template <typename T>
constexpr Vec3<T> min(const Vec3<T>& a, const Vec3<T>& b) {
  return {std::min(a.x, b.x), std::min(a.y, b.y), std::min(a.z, b.z)};
}

template <typename T>
constexpr Vec3<T> max(const Vec3<T>& a, const Vec3<T>& b) {
  return {std::max(a.x, b.x), std::max(a.y, b.y), std::max(a.z, b.z)};
}

template <typename T>
Vec3<T> abs(const Vec3<T>& a) {
  return {std::abs(a.x), std::abs(a.y), std::abs(a.z)};
}

/// Index of the component with the largest value.
template <typename T>
constexpr int max_dimension(const Vec3<T>& a) {
  if (a.x >= a.y && a.x >= a.z) {
    return 0;
  }
  return a.y >= a.z ? 1 : 2;
}

/// Mirror `d` about the plane with unit normal `n`.
template <typename T>
constexpr Vec3<T> reflect(const Vec3<T>& d, const Vec3<T>& n) {
  return d - n * (T(2) * dot(d, n));
}

}  // namespace bsar
