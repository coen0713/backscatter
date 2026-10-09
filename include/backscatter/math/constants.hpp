#pragma once

#include <numbers>

namespace bsar {

inline constexpr double kPi = std::numbers::pi;
inline constexpr double kSpeedOfLight = 299'792'458.0;  // m/s, exact

constexpr double deg_to_rad(double deg) { return deg * kPi / 180.0; }

constexpr double rad_to_deg(double rad) { return rad * 180.0 / kPi; }

}  // namespace bsar
