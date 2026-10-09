#pragma once

#include <utility>

#include "backscatter/math/aabb.hpp"
#include "backscatter/math/vec3.hpp"
#include "backscatter/sensor/geodesy.hpp"
#include "backscatter/sensor/orbit.hpp"

namespace bsar {

/// Sensor antenna phase centre as a function of time, in scene-local ENU
/// coordinates (metres, double precision).
class Trajectory {
 public:
  Trajectory() = default;
  Trajectory(const Trajectory&) = default;
  Trajectory& operator=(const Trajectory&) = default;
  Trajectory(Trajectory&&) = default;
  Trajectory& operator=(Trajectory&&) = default;
  virtual ~Trajectory() = default;

  [[nodiscard]] virtual Vec3d position(double t) const = 0;
  [[nodiscard]] virtual Vec3d velocity(double t) const = 0;
  [[nodiscard]] virtual double t_begin() const = 0;
  [[nodiscard]] virtual double t_end() const = 0;
};

/// Constant-velocity straight line: p(t) = p0 + v t.
class LinearTrajectory final : public Trajectory {
 public:
  LinearTrajectory(const Vec3d& p0, const Vec3d& v, double t_begin, double t_end)
      : p0_(p0), v_(v), t0_(t_begin), t1_(t_end) {}

  [[nodiscard]] Vec3d position(double t) const override { return p0_ + v_ * t; }
  [[nodiscard]] Vec3d velocity(double /*t*/) const override { return v_; }
  [[nodiscard]] double t_begin() const override { return t0_; }
  [[nodiscard]] double t_end() const override { return t1_; }

 private:
  Vec3d p0_;
  Vec3d v_;
  double t0_;
  double t1_;
};

/// An Earth-fixed orbit expressed in a scene's ENU frame. Earth-fixed
/// velocities make zero-Doppler geometry correct without extra terms for
/// Earth rotation.
class OrbitTrajectory final : public Trajectory {
 public:
  OrbitTrajectory(Orbit orbit, const EnuFrame& frame) : orbit_(std::move(orbit)), frame_(frame) {}

  [[nodiscard]] Vec3d position(double t) const override {
    return frame_.to_enu(orbit_.position(t));
  }
  [[nodiscard]] Vec3d velocity(double t) const override {
    return frame_.rotate_to_enu(orbit_.velocity(t));
  }
  [[nodiscard]] double t_begin() const override { return orbit_.t_begin(); }
  [[nodiscard]] double t_end() const override { return orbit_.t_end(); }

 private:
  Orbit orbit_;
  EnuFrame frame_;
};

enum class LookSide { Right, Left };

/// Flat-earth side-looking geometry: a straight track heading north (+y) at
/// `altitude`, offset east/west so the scene centre is seen at
/// `incidence_deg`. Passes the scene centre at t = 0.
LinearTrajectory make_side_looking_track(const Aabb& scene, double altitude, double speed,
                                         double incidence_deg, LookSide side);

}  // namespace bsar
