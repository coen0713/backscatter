#include "backscatter/sensor/orbit.hpp"

#include <algorithm>
#include <charconv>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace bsar {
namespace {

// Howard Hinnant's days_from_civil: days since 1970-01-01 in the proleptic
// Gregorian calendar.
long long days_from_civil(long long y, unsigned m, unsigned d) {
  y -= m <= 2 ? 1 : 0;
  const long long era = (y >= 0 ? y : y - 399) / 400;
  const auto yoe = static_cast<unsigned>(y - era * 400);
  const unsigned doy = (153 * (m + (m > 2 ? -3 : 9)) + 2) / 5 + d - 1;
  const unsigned doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
  return era * 146097 + static_cast<long long>(doe) - 719468;
}

double parse_double(std::string_view s) {
  while (!s.empty() && (s.front() == ' ' || s.front() == '\n' || s.front() == '\r' ||
                        s.front() == '\t' || s.front() == '+')) {
    s.remove_prefix(1);
  }
  // std::from_chars for double is not available on every standard library we
  // target, so go through strtod on a bounded copy.
  const std::string copy(s);
  char* end = nullptr;
  const double v = std::strtod(copy.c_str(), &end);
  if (end == copy.c_str()) {
    throw std::runtime_error("orbit: cannot parse number '" + copy + "'");
  }
  return v;
}

int parse_int(std::string_view s) {
  int v = 0;
  const auto [ptr, ec] = std::from_chars(s.data(), s.data() + s.size(), v);
  if (ec != std::errc{} || ptr != s.data() + s.size()) {
    throw std::runtime_error("orbit: cannot parse integer '" + std::string(s) + "'");
  }
  return v;
}

/// Content of the first <tag ...>content</tag> inside `block`.
std::string_view extract_tag(std::string_view block, std::string_view tag) {
  std::size_t pos = 0;
  while (true) {
    pos = block.find('<', pos);
    if (pos == std::string_view::npos) {
      throw std::runtime_error("orbit: missing <" + std::string(tag) + ">");
    }
    const std::string_view rest = block.substr(pos + 1);
    if (rest.substr(0, tag.size()) == tag && rest.size() > tag.size() &&
        (rest[tag.size()] == '>' || rest[tag.size()] == ' ')) {
      break;
    }
    ++pos;
  }
  const std::size_t open_end = block.find('>', pos);
  const std::string close = "</" + std::string(tag) + ">";
  const std::size_t close_pos = block.find(close, open_end);
  if (open_end == std::string_view::npos || close_pos == std::string_view::npos) {
    throw std::runtime_error("orbit: malformed <" + std::string(tag) + ">");
  }
  return block.substr(open_end + 1, close_pos - open_end - 1);
}

}  // namespace

double parse_utc(std::string_view text) {
  while (!text.empty() && (text.front() == ' ' || text.front() == '\n')) {
    text.remove_prefix(1);
  }
  if (const auto eq = text.find('='); eq != std::string_view::npos && eq < 4) {
    text.remove_prefix(eq + 1);
  }
  if (text.size() < 19 || text[4] != '-' || text[7] != '-' || text[10] != 'T' || text[13] != ':' ||
      text[16] != ':') {
    throw std::runtime_error("parse_utc: expected YYYY-MM-DDThh:mm:ss, got '" + std::string(text) +
                             "'");
  }
  const int year = parse_int(text.substr(0, 4));
  const int month = parse_int(text.substr(5, 2));
  const int day = parse_int(text.substr(8, 2));
  const int hour = parse_int(text.substr(11, 2));
  const int minute = parse_int(text.substr(14, 2));
  std::string_view sec_text = text.substr(17);
  const auto end = sec_text.find_first_not_of("0123456789.");
  if (end != std::string_view::npos) {
    sec_text = sec_text.substr(0, end);
  }
  const double second = parse_double(sec_text);

  const long long days =
      days_from_civil(year, static_cast<unsigned>(month), static_cast<unsigned>(day)) -
      days_from_civil(2000, 1, 1);
  return static_cast<double>(days) * 86400.0 + hour * 3600.0 + minute * 60.0 + second;
}

Orbit::Orbit(std::vector<StateVector> vectors) : vectors_(std::move(vectors)) {
  std::sort(vectors_.begin(), vectors_.end(),
            [](const StateVector& a, const StateVector& b) { return a.t < b.t; });
  vectors_.erase(std::unique(vectors_.begin(), vectors_.end(),
                             [](const StateVector& a, const StateVector& b) { return a.t == b.t; }),
                 vectors_.end());
  if (vectors_.size() < 2) {
    throw std::invalid_argument("Orbit: need at least two distinct state vectors");
  }
}

double Orbit::t_begin() const { return vectors_.front().t; }

double Orbit::t_end() const { return vectors_.back().t; }

std::size_t Orbit::segment(double t) const {
  const auto it = std::upper_bound(vectors_.begin(), vectors_.end(), t,
                                   [](double v, const StateVector& s) { return v < s.t; });
  const auto idx = static_cast<std::size_t>(std::max<std::ptrdiff_t>(it - vectors_.begin(), 1));
  return std::min(idx, vectors_.size() - 1) - 1;
}

Vec3d Orbit::position(double t) const {
  const std::size_t i = segment(t);
  const StateVector& a = vectors_[i];
  const StateVector& b = vectors_[i + 1];
  const double h = b.t - a.t;
  const double s = (t - a.t) / h;
  const double s2 = s * s;
  const double s3 = s2 * s;
  const double h00 = 2 * s3 - 3 * s2 + 1;
  const double h10 = s3 - 2 * s2 + s;
  const double h01 = -2 * s3 + 3 * s2;
  const double h11 = s3 - s2;
  return a.position * h00 + a.velocity * (h10 * h) + b.position * h01 + b.velocity * (h11 * h);
}

Vec3d Orbit::velocity(double t) const {
  const std::size_t i = segment(t);
  const StateVector& a = vectors_[i];
  const StateVector& b = vectors_[i + 1];
  const double h = b.t - a.t;
  const double s = (t - a.t) / h;
  const double s2 = s * s;
  const double d00 = 6 * s2 - 6 * s;
  const double d10 = 3 * s2 - 4 * s + 1;
  const double d01 = -6 * s2 + 6 * s;
  const double d11 = 3 * s2 - 2 * s;
  return (a.position * d00 + b.position * d01) / h + a.velocity * d10 + b.velocity * d11;
}

Orbit parse_eof(std::string_view xml) {
  std::vector<StateVector> vectors;
  std::size_t pos = 0;
  while ((pos = xml.find("<OSV>", pos)) != std::string_view::npos) {
    const std::size_t end = xml.find("</OSV>", pos);
    if (end == std::string_view::npos) {
      throw std::runtime_error("parse_eof: unterminated <OSV>");
    }
    const std::string_view block = xml.substr(pos, end - pos);
    StateVector sv;
    sv.t = parse_utc(extract_tag(block, "UTC"));
    sv.position = {parse_double(extract_tag(block, "X")), parse_double(extract_tag(block, "Y")),
                   parse_double(extract_tag(block, "Z"))};
    sv.velocity = {parse_double(extract_tag(block, "VX")), parse_double(extract_tag(block, "VY")),
                   parse_double(extract_tag(block, "VZ"))};
    vectors.push_back(sv);
    pos = end;
  }
  return Orbit(std::move(vectors));
}

Orbit load_eof(const std::filesystem::path& path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) {
    throw std::runtime_error("load_eof: cannot open " + path.string());
  }
  std::ostringstream ss;
  ss << in.rdbuf();
  return parse_eof(ss.str());
}

}  // namespace bsar
