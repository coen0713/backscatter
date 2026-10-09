#include "backscatter/trace/geometric.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

#include "backscatter/sensor/sar_geometry.hpp"
#include "backscatter/trace/sensor_rays.hpp"
#include "backscatter/util/parallel.hpp"

namespace bsar {

double find_zero_doppler_time(const Trajectory& trajectory, const Vec3d& target) {
  constexpr int kSamples = 4096;
  const double t0 = trajectory.t_begin();
  const double t1 = trajectory.t_end();
  double best_t = 0.5 * (t0 + t1);
  double best_r = std::numeric_limits<double>::infinity();
  for (int i = 0; i <= kSamples; ++i) {
    const double t = t0 + (t1 - t0) * i / kSamples;
    const double r = length(trajectory.position(t) - target);
    if (r < best_r) {
      best_r = r;
      best_t = t;
    }
  }
  return solve_zero_doppler(trajectory, target, best_t).t;
}

namespace {

struct FirstHit {
  bool valid = false;
  double range = 0.0;  // one-way slant range [m]
  Vec3d position;
  double cos_incidence = 0.0;
  double ray_range = 0.0;  // distance along the ray, for the spacing estimate
};

/// Spread `value` uniformly over [center - width/2, center + width/2] (in
/// bin units) and add it to both rows.
void splat(std::vector<double>& row, std::vector<double>& bounce_row, double center, double width,
           double value) {
  const auto bins = static_cast<double>(row.size());
  if (width < 1e-6) {
    const double k = std::floor(center);
    if (k >= 0.0 && k < bins) {
      row[static_cast<std::size_t>(k)] += value;
      bounce_row[static_cast<std::size_t>(k)] += value;
    }
    return;
  }
  const double a = center - 0.5 * width;
  const double b = center + 0.5 * width;
  const double density = value / width;
  for (double k = std::max(0.0, std::floor(a)); k < std::min(bins, b); k += 1.0) {
    const double overlap = std::min(b, k + 1.0) - std::max(a, k);
    if (overlap > 0.0) {
      const auto i = static_cast<std::size_t>(k);
      row[i] += density * overlap;
      bounce_row[i] += density * overlap;
    }
  }
}

/// Add `len` of coverage multiplicity `mult` over [a, b) to the per-bin
/// accumulators (b0 = bins with no coverage, b2 = bins covered twice or more).
void add_measure(double a, double b, int mult, double near, double dr, std::vector<double>& m0,
                 std::vector<double>& m2) {
  if (b <= a || (mult != 0 && mult < 2)) {
    return;
  }
  const auto n = static_cast<std::ptrdiff_t>(m0.size());
  auto first = static_cast<std::ptrdiff_t>(std::floor((a - near) / dr));
  auto last = static_cast<std::ptrdiff_t>(std::floor((b - near) / dr));
  first = std::max<std::ptrdiff_t>(first, 0);
  last = std::min<std::ptrdiff_t>(last, n - 1);
  for (std::ptrdiff_t k = first; k <= last; ++k) {
    const double lo = near + static_cast<double>(k) * dr;
    const double overlap = std::min(b, lo + dr) - std::max(a, lo);
    if (overlap <= 0.0) {
      continue;
    }
    (mult == 0 ? m0 : m2)[static_cast<std::size_t>(k)] += overlap;
  }
}

/// Layover and shadow for one azimuth line from the sequence of first hits
/// (ordered by look angle). Each pair of neighbouring hits on a continuous
/// surface covers the slant-range interval between them; bins covered more
/// than once are in layover, bins inside the swath covered by nothing are in
/// radar shadow.
void compute_masks(const std::vector<FirstHit>& hits, double dtheta, double jump_factor,
                   double near, double dr, std::uint8_t* layover, std::uint8_t* shadow,
                   std::size_t bins) {
  std::vector<std::pair<double, int>> events;
  double span_lo = std::numeric_limits<double>::infinity();
  double span_hi = -std::numeric_limits<double>::infinity();
  for (std::size_t i = 0; i < hits.size(); ++i) {
    if (!hits[i].valid) {
      continue;
    }
    span_lo = std::min(span_lo, hits[i].range);
    span_hi = std::max(span_hi, hits[i].range);
    if (i + 1 >= hits.size() || !hits[i + 1].valid) {
      continue;
    }
    const FirstHit& a = hits[i];
    const FirstHit& b = hits[i + 1];
    const double c = std::max({a.cos_incidence, b.cos_incidence, 0.05});
    const double expected = a.ray_range * dtheta / c;
    if (length(b.position - a.position) > jump_factor * expected) {
      continue;  // occlusion boundary: nothing between these hits is seen
    }
    const double lo = std::min(a.range, b.range);
    const double hi = std::max(a.range, b.range);
    if (hi > lo) {
      events.emplace_back(lo, +1);
      events.emplace_back(hi, -1);
    }
  }
  if (!(span_hi > span_lo)) {
    return;
  }
  std::sort(events.begin(), events.end());
  std::vector<double> m0(bins, 0.0);
  std::vector<double> m2(bins, 0.0);
  double prev = span_lo;
  int mult = 0;
  for (const auto& [x, delta] : events) {
    add_measure(prev, x, mult, near, dr, m0, m2);
    mult += delta;
    prev = x;
  }
  add_measure(prev, span_hi, mult, near, dr, m0, m2);
  for (std::size_t k = 0; k < bins; ++k) {
    layover[k] = m2[k] > 0.05 * dr ? 1 : 0;
    shadow[k] = m0[k] > 0.5 * dr ? 1 : 0;
  }
}

}  // namespace

GeometricImage render_geometric(const Scene& scene, const Trajectory& trajectory,
                                const GeometricConfig& config) {
  const FacetScatteringModel model;
  return render_geometric(scene, trajectory, config, model);
}

GeometricImage render_geometric(const Scene& scene, const Trajectory& trajectory,
                                const GeometricConfig& config, const ScatteringModel& model) {
  if (!scene.built()) {
    throw std::invalid_argument("render_geometric: scene is not built");
  }
  if (config.range_spacing <= 0.0 || config.azimuth_spacing <= 0.0 || config.rays_per_bin <= 0.0 ||
      config.max_bounces < 1) {
    throw std::invalid_argument("render_geometric: invalid configuration");
  }
  const Aabb& bounds = scene.bounds();
  const double scene_extent = length(Vec3d(bounds.extent()));
  const double eps = 1e-3 + 1e-6 * scene_extent;  // ray offset off surfaces [m]

  // --- Acquisition window from the scene corners -------------------------
  double t_lo = std::numeric_limits<double>::infinity();
  double t_hi = -t_lo;
  double r_lo = std::numeric_limits<double>::infinity();
  double r_hi = -r_lo;
  for (int c = 0; c < 8; ++c) {
    const Vec3d corner(bounds.corner(c));
    const double t = find_zero_doppler_time(trajectory, corner);
    t_lo = std::min(t_lo, t);
    t_hi = std::max(t_hi, t);
    const double r = length(trajectory.position(t) - corner);
    r_lo = std::min(r_lo, r);
    r_hi = std::max(r_hi, r);
  }
  const double speed = length(trajectory.velocity(0.5 * (t_lo + t_hi)));
  const double line_interval = config.azimuth_spacing / speed;
  const double t_start = std::isnan(config.t_start) ? t_lo - line_interval : config.t_start;
  const double t_end = std::isnan(config.t_end) ? t_hi + line_interval : config.t_end;
  const double dr = config.range_spacing;
  const double near = std::isnan(config.near_range) ? r_lo - 2 * dr : config.near_range;
  const double far = std::isnan(config.far_range) ? r_hi + 2 * dr : config.far_range;
  if (!(t_end > t_start) || !(far > near)) {
    throw std::invalid_argument("render_geometric: empty acquisition window");
  }

  const auto lines = static_cast<std::size_t>(std::floor((t_end - t_start) / line_interval)) + 1;
  const auto bins = static_cast<std::size_t>(std::ceil((far - near) / dr));
  const auto max_bounces = static_cast<std::size_t>(config.max_bounces);

  GeometricImage out;
  out.intensity = Image<float>(bins, lines);
  out.bounce.assign(max_bounces, Image<float>(bins, lines));
  out.layover = Image<std::uint8_t>(bins, lines);
  out.shadow = Image<std::uint8_t>(bins, lines);
  out.near_range = near;
  out.range_spacing = dr;
  out.t_start = t_start;
  out.line_interval = line_interval;

  std::vector<std::size_t> rays_per_line(lines, 0);
  const Vec3d center(bounds.center());
  const double norm = 1.0 / (dr * config.azimuth_spacing);

  parallel_for(0, lines, config.threads, [&](std::size_t line) {
    const double t = t_start + static_cast<double>(line) * line_interval;
    const Vec3d p = trajectory.position(t);
    const Vec3d ev = normalize(trajectory.velocity(t));
    // Zero-Doppler plane basis: e1 towards the scene centre, e2 = ev x e1.
    Vec3d to_center = center - p;
    to_center -= ev * dot(to_center, ev);
    const Vec3d e1 = normalize(to_center);
    const Vec3d e2 = cross(ev, e1);

    double th_lo = std::numeric_limits<double>::infinity();
    double th_hi = -th_lo;
    for (int c = 0; c < 8; ++c) {
      const Vec3d d = Vec3d(bounds.corner(c)) - p;
      const double th = std::atan2(dot(d, e2), dot(d, e1));
      th_lo = std::min(th_lo, th);
      th_hi = std::max(th_hi, th);
    }
    const double r_mid = length(to_center);
    const double dtheta = dr / (config.rays_per_bin * r_mid);
    const auto n_rays = static_cast<std::size_t>(std::ceil((th_hi - th_lo) / dtheta)) + 1;

    std::vector<double> row(bins, 0.0);
    std::vector<std::vector<double>> bounce_rows(max_bounces, std::vector<double>(bins, 0.0));
    std::vector<FirstHit> first_hits(n_rays);

    for (std::size_t i = 0; i < n_rays; ++i) {
      const double th = th_lo + (static_cast<double>(i) + 0.5) * dtheta;
      Vec3d dir = e1 * std::cos(th) + e2 * std::sin(th);
      const auto clipped = clip_to_bounds(p, dir, bounds, 1.0);
      if (!clipped) {
        continue;
      }
      Ray ray = clipped->ray;
      double path = clipped->t_offset;
      double weight = 1.0;
      double tube_area = 0.0;  // cross-section of this ray tube [m^2]

      for (std::size_t b = 0; b < max_bounces; ++b) {
        Hit hit;
        if (!scene.intersect(ray, hit)) {
          break;
        }
        path += hit.t;
        const Vec3d pos = Vec3d(ray.origin) + Vec3d(ray.dir) * static_cast<double>(hit.t);
        Vec3d n(scene.mesh.normal(hit.prim));
        if (dot(n, dir) > 0.0) {
          n = -n;  // surfaces are two-sided
        }
        const Material& mat = scene.material_of(hit.prim);
        const double cos_in = std::max(-dot(n, dir), 1e-6);
        const Vec3d to_sensor = p - pos;
        const double dist_sensor = length(to_sensor);
        const Vec3d s = to_sensor / dist_sensor;

        double sigma0 = 0.0;
        if (b == 0) {
          tube_area = path * dtheta * config.azimuth_spacing;
          sigma0 = model.backscatter(mat, cos_in, config.polarization);
          first_hits[i] = {true, path, pos, cos_in, path};
        } else if (dot(n, s) > 0.0) {
          const Ray shadow_ray{Vec3f(pos + n * eps), Vec3f(s), 0.0f,
                               std::numeric_limits<float>::infinity()};
          if (!scene.occluded(shadow_ray)) {
            sigma0 = model.bistatic(mat, n, -dir, s, config.polarization);
          }
        }
        if (sigma0 > 0.0) {
          const double contrib = weight * sigma0 * tube_area / cos_in * norm;
          const double r_eq = 0.5 * (path + dist_sensor);
          // A direct return's footprint spans R dtheta tan(theta_local) in
          // slant range; spreading it over that extent (instead of point
          // binning) removes aliasing between ray and bin spacings.
          // Multi-bounce paths stay point-binned: for corner reflectors every
          // path has the same length, so they focus to one range.
          double footprint = 0.0;
          if (b == 0) {
            const double sin_in = std::sqrt(std::max(0.0, 1.0 - cos_in * cos_in));
            footprint = std::min(path * dtheta * sin_in / cos_in, 4.0 * dr);
          }
          splat(row, bounce_rows[b], (r_eq - near) / dr, footprint / dr, contrib);
        }

        weight *= model.specular_reflectance(mat, cos_in, config.polarization);
        if (weight < 1e-6) {
          break;
        }
        dir = reflect(dir, n);
        ray = Ray{Vec3f(pos + n * eps), Vec3f(dir), 0.0f, std::numeric_limits<float>::infinity()};
      }
    }

    float* out_row = out.intensity.row(line);
    for (std::size_t k = 0; k < bins; ++k) {
      out_row[k] = static_cast<float>(row[k]);
    }
    for (std::size_t b = 0; b < max_bounces; ++b) {
      float* br = out.bounce[b].row(line);
      for (std::size_t k = 0; k < bins; ++k) {
        br[k] = static_cast<float>(bounce_rows[b][k]);
      }
    }
    compute_masks(first_hits, dtheta, config.occlusion_jump_factor, near, dr, out.layover.row(line),
                  out.shadow.row(line), bins);
    rays_per_line[line] = n_rays;
  });

  for (const std::size_t n : rays_per_line) {
    out.rays_traced += n;
  }
  return out;
}

}  // namespace bsar
