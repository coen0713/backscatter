#include "backscatter/scene/material.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace bsar {

namespace materials {

// Permittivities: Hallikainen et al. (1985) soil model at 5 GHz; Debye model
// for fresh water at 20 C; Ulaby & Long (2014) tables for building materials.
// Diffuse coefficients are chosen so sigma0 at 35 deg incidence falls in the
// typical C-band ranges of the same references.

Material dry_soil() { return {"dry_soil", {3.5, -0.2}, 0.08, 2.0, 0.15, false}; }

Material wet_soil() { return {"wet_soil", {15.0, -3.0}, 0.15, 2.0, 0.15, false}; }

Material water() { return {"water", {68.0, -21.0}, 0.003, 2.0, 0.05, false}; }

Material concrete() { return {"concrete", {5.5, -0.6}, 0.06, 1.5, 0.10, false}; }

Material metal() { return {"metal", {1.0, 0.0}, 0.01, 1.5, 0.05, true}; }

Material asphalt() { return {"asphalt", {4.0, -0.3}, 0.02, 2.0, 0.08, false}; }

Material by_name(const std::string& name) {
  for (const auto& m : {dry_soil(), wet_soil(), water(), concrete(), metal(), asphalt()}) {
    if (m.name == name) {
      return m;
    }
  }
  throw std::invalid_argument("unknown material '" + name + "'");
}

}  // namespace materials

std::complex<double> fresnel_reflection(const Material& m, double cos_theta, Polarization pol) {
  if (m.perfect_conductor) {
    return pol == Polarization::HH ? -1.0 : 1.0;
  }
  const double c = std::clamp(cos_theta, 0.0, 1.0);
  const double s2 = 1.0 - c * c;
  const std::complex<double> eps = m.permittivity;
  const std::complex<double> root = std::sqrt(eps - s2);
  if (pol == Polarization::HH) {
    return (c - root) / (c + root);
  }
  return (eps * c - root) / (eps * c + root);
}

namespace {

double geometric_optics_lobe(double cos_beta, double rms_slope) {
  const double c = std::max(cos_beta, 1e-3);
  const double c2 = c * c;
  const double tan2 = (1.0 - c2) / c2;
  const double s2 = rms_slope * rms_slope;
  return std::exp(-tan2 / (2.0 * s2)) / (2.0 * s2 * c2 * c2);
}

}  // namespace

double FacetScatteringModel::backscatter(const Material& m, double cos_theta,
                                         Polarization pol) const {
  if (cos_theta <= 0.0) {
    return 0.0;
  }
  const double diffuse = m.diffuse_coefficient * std::pow(cos_theta, m.diffuse_exponent);
  const double r0 = std::norm(fresnel_reflection(m, 1.0, pol));
  return diffuse + r0 * geometric_optics_lobe(cos_theta, m.rms_slope);
}

double FacetScatteringModel::bistatic(const Material& m, const Vec3d& n, const Vec3d& to_source,
                                      const Vec3d& to_receiver, Polarization pol) const {
  const double ci = dot(n, to_source);
  const double cs = dot(n, to_receiver);
  if (ci <= 0.0 || cs <= 0.0) {
    return 0.0;
  }
  const double diffuse = m.diffuse_coefficient * std::pow(ci * cs, 0.5 * m.diffuse_exponent);
  const Vec3d h = normalize(to_source + to_receiver);
  const double r = std::norm(fresnel_reflection(m, dot(to_source, h), pol));
  return diffuse + r * geometric_optics_lobe(dot(n, h), m.rms_slope);
}

double FacetScatteringModel::specular_reflectance(const Material& m, double cos_theta,
                                                  Polarization pol) const {
  return std::norm(fresnel_reflection(m, cos_theta, pol));
}

}  // namespace bsar
