#pragma once

#include <complex>
#include <string>

#include "backscatter/math/vec3.hpp"

namespace bsar {

enum class Polarization { HH, VV };

/// Surface material for microwave scattering. Values are C-band (~5.4 GHz)
/// order-of-magnitude figures; see docs/scattering.md for sources.
struct Material {
  std::string name;
  std::complex<double> permittivity{1.0, 0.0};  // relative, eps' - j eps''
  double diffuse_coefficient = 0.0;             // k in sigma0_d = k cos^n(theta)
  double diffuse_exponent = 2.0;                // n
  double rms_slope = 0.1;                       // s, geometric-optics facet slope
  bool perfect_conductor = false;
};

namespace materials {
Material dry_soil();
Material wet_soil();
Material water();
Material concrete();
Material metal();
Material asphalt();
/// Look up one of the presets above by name; throws on unknown names.
Material by_name(const std::string& name);
}  // namespace materials

/// Fresnel amplitude reflection coefficient for a planar interface between
/// vacuum and a medium of complex relative permittivity, at incidence angle
/// acos(cos_theta). A perfect conductor returns -1 (H) / +1 (V).
std::complex<double> fresnel_reflection(const Material& m, double cos_theta, Polarization pol);

/// The single interface between integrators and surface physics. A better
/// model (e.g. a simplified IEM) can replace the default without touching the
/// integrators.
class ScatteringModel {
 public:
  ScatteringModel() = default;
  ScatteringModel(const ScatteringModel&) = default;
  ScatteringModel& operator=(const ScatteringModel&) = default;
  ScatteringModel(ScatteringModel&&) = default;
  ScatteringModel& operator=(ScatteringModel&&) = default;
  virtual ~ScatteringModel() = default;

  /// Monostatic backscatter coefficient sigma0 (m^2/m^2) at local incidence
  /// acos(cos_theta).
  [[nodiscard]] virtual double backscatter(const Material& m, double cos_theta,
                                           Polarization pol) const = 0;

  /// Bistatic scattering coefficient for light arriving along `to_source`
  /// (unit, pointing away from the surface) and leaving along `to_receiver`,
  /// with unit surface normal `n`. Reduces to backscatter() when the two
  /// directions coincide.
  [[nodiscard]] virtual double bistatic(const Material& m, const Vec3d& n, const Vec3d& to_source,
                                        const Vec3d& to_receiver, Polarization pol) const = 0;

  /// Fraction of power carried by the specular (mirror) reflection.
  [[nodiscard]] virtual double specular_reflectance(const Material& m, double cos_theta,
                                                    Polarization pol) const = 0;
};

/// Default model: a diffuse term sigma0_d = k cos^n(theta) plus a
/// Kirchhoff geometric-optics (stationary phase) specular lobe
///   sigma0_s = |R(0)|^2 exp(-tan^2(beta) / (2 s^2)) / (2 s^2 cos^4(beta)),
/// where beta is the angle between the normal and the bistatic half vector
/// (Ulaby & Long, "Microwave Radar and Radiometric Remote Sensing", 2014).
class FacetScatteringModel final : public ScatteringModel {
 public:
  [[nodiscard]] double backscatter(const Material& m, double cos_theta,
                                   Polarization pol) const override;
  [[nodiscard]] double bistatic(const Material& m, const Vec3d& n, const Vec3d& to_source,
                                const Vec3d& to_receiver, Polarization pol) const override;
  [[nodiscard]] double specular_reflectance(const Material& m, double cos_theta,
                                            Polarization pol) const override;
};

}  // namespace bsar
