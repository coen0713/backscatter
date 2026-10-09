#include "backscatter/image/io.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <fstream>
#include <stdexcept>

namespace bsar {
namespace {

std::array<std::uint32_t, 256> make_crc_table() {
  std::array<std::uint32_t, 256> table{};
  for (std::uint32_t n = 0; n < 256; ++n) {
    std::uint32_t c = n;
    for (int k = 0; k < 8; ++k) {
      c = (c & 1u) != 0 ? 0xEDB88320u ^ (c >> 1) : c >> 1;
    }
    table[n] = c;
  }
  return table;
}

std::uint32_t crc32(const std::uint8_t* data, std::size_t n, std::uint32_t crc = 0) {
  static const auto table = make_crc_table();
  crc = ~crc;
  for (std::size_t i = 0; i < n; ++i) {
    crc = table[(crc ^ data[i]) & 0xFFu] ^ (crc >> 8);
  }
  return ~crc;
}

void put_u32_be(std::vector<std::uint8_t>& out, std::uint32_t v) {
  out.push_back(static_cast<std::uint8_t>(v >> 24));
  out.push_back(static_cast<std::uint8_t>(v >> 16));
  out.push_back(static_cast<std::uint8_t>(v >> 8));
  out.push_back(static_cast<std::uint8_t>(v));
}

void write_chunk(std::ofstream& f, const char* type, const std::vector<std::uint8_t>& data) {
  std::vector<std::uint8_t> buf;
  put_u32_be(buf, static_cast<std::uint32_t>(data.size()));
  buf.insert(buf.end(), type, type + 4);
  buf.insert(buf.end(), data.begin(), data.end());
  const std::uint32_t crc = crc32(buf.data() + 4, buf.size() - 4);
  put_u32_be(buf, crc);
  f.write(reinterpret_cast<const char*>(buf.data()), static_cast<std::streamsize>(buf.size()));
}

std::ofstream open_binary(const std::filesystem::path& path) {
  std::ofstream f(path, std::ios::binary);
  if (!f) {
    throw std::runtime_error("cannot write " + path.string());
  }
  return f;
}

void write_npy_raw(const std::filesystem::path& path, const char* descr, std::size_t height,
                   std::size_t width, const void* data, std::size_t bytes) {
  std::string header = "{'descr': '" + std::string(descr) +
                       "', 'fortran_order': False, 'shape': (" + std::to_string(height) + ", " +
                       std::to_string(width) + "), }";
  const std::size_t total = 10 + header.size() + 1;
  header.append((64 - total % 64) % 64, ' ');
  header.push_back('\n');
  auto f = open_binary(path);
  const char magic[] = "\x93NUMPY\x01\x00";
  f.write(magic, 8);
  const auto len = static_cast<std::uint16_t>(header.size());
  const std::array<char, 2> len_bytes{static_cast<char>(len & 0xFF), static_cast<char>(len >> 8)};
  f.write(len_bytes.data(), 2);
  f.write(header.data(), static_cast<std::streamsize>(header.size()));
  f.write(static_cast<const char*>(data), static_cast<std::streamsize>(bytes));
}

}  // namespace

void write_png(const std::filesystem::path& path, std::size_t width, std::size_t height,
               int channels, const std::vector<std::uint8_t>& pixels) {
  if (channels != 1 && channels != 3) {
    throw std::invalid_argument("write_png: channels must be 1 or 3");
  }
  const std::size_t stride = width * static_cast<std::size_t>(channels);
  if (pixels.size() != stride * height) {
    throw std::invalid_argument("write_png: pixel buffer has the wrong size");
  }
  // Raw scanlines, each prefixed with filter type 0.
  std::vector<std::uint8_t> raw;
  raw.reserve((stride + 1) * height);
  for (std::size_t y = 0; y < height; ++y) {
    raw.push_back(0);
    raw.insert(raw.end(), pixels.begin() + static_cast<std::ptrdiff_t>(y * stride),
               pixels.begin() + static_cast<std::ptrdiff_t>((y + 1) * stride));
  }
  // zlib stream of stored blocks.
  std::vector<std::uint8_t> z{0x78, 0x01};
  std::size_t pos = 0;
  do {
    const std::size_t n = std::min<std::size_t>(65535, raw.size() - pos);
    const bool last = pos + n == raw.size();
    z.push_back(last ? 1 : 0);
    z.push_back(static_cast<std::uint8_t>(n & 0xFF));
    z.push_back(static_cast<std::uint8_t>(n >> 8));
    z.push_back(static_cast<std::uint8_t>(~n & 0xFF));
    z.push_back(static_cast<std::uint8_t>((~n >> 8) & 0xFF));
    z.insert(z.end(), raw.begin() + static_cast<std::ptrdiff_t>(pos),
             raw.begin() + static_cast<std::ptrdiff_t>(pos + n));
    pos += n;
  } while (pos < raw.size());
  std::uint32_t a = 1;
  std::uint32_t b = 0;
  for (const std::uint8_t byte : raw) {
    a = (a + byte) % 65521u;
    b = (b + a) % 65521u;
  }
  put_u32_be(z, (b << 16) | a);

  std::vector<std::uint8_t> ihdr;
  put_u32_be(ihdr, static_cast<std::uint32_t>(width));
  put_u32_be(ihdr, static_cast<std::uint32_t>(height));
  ihdr.insert(ihdr.end(), {8, static_cast<std::uint8_t>(channels == 1 ? 0 : 2), 0, 0, 0});

  auto f = open_binary(path);
  const std::array<char, 8> signature{'\x89', 'P', 'N', 'G', '\r', '\n', '\x1a', '\n'};
  f.write(signature.data(), 8);
  write_chunk(f, "IHDR", ihdr);
  write_chunk(f, "IDAT", z);
  write_chunk(f, "IEND", {});
}

void write_png(const std::filesystem::path& path, const Image<std::uint8_t>& gray) {
  write_png(path, gray.width, gray.height, 1, gray.data);
}

void write_npy(const std::filesystem::path& path, const Image<float>& image) {
  write_npy_raw(path, "<f4", image.height, image.width, image.data.data(),
                image.data.size() * sizeof(float));
}

void write_npy(const std::filesystem::path& path, const Image<std::complex<float>>& image) {
  write_npy_raw(path, "<c8", image.height, image.width, image.data.data(),
                image.data.size() * sizeof(std::complex<float>));
}

void write_npy(const std::filesystem::path& path, const Image<std::uint8_t>& image) {
  write_npy_raw(path, "|u1", image.height, image.width, image.data.data(), image.data.size());
}

NpyArray read_npy(const std::filesystem::path& path) {
  std::ifstream f(path, std::ios::binary);
  if (!f) {
    throw std::runtime_error("read_npy: cannot open " + path.string());
  }
  std::array<char, 8> magic{};
  f.read(magic.data(), 8);
  if (std::memcmp(magic.data(), "\x93NUMPY", 6) != 0) {
    throw std::runtime_error("read_npy: not a .npy file: " + path.string());
  }
  std::size_t header_len = 0;
  if (magic[6] == 1) {
    std::array<unsigned char, 2> b{};
    f.read(reinterpret_cast<char*>(b.data()), 2);
    header_len = b[0] | (static_cast<std::size_t>(b[1]) << 8);
  } else {
    std::array<unsigned char, 4> b{};
    f.read(reinterpret_cast<char*>(b.data()), 4);
    header_len = b[0] | (static_cast<std::size_t>(b[1]) << 8) |
                 (static_cast<std::size_t>(b[2]) << 16) | (static_cast<std::size_t>(b[3]) << 24);
  }
  std::string header(header_len, '\0');
  f.read(header.data(), static_cast<std::streamsize>(header_len));

  NpyArray arr;
  const auto d = header.find("'descr':");
  const auto q1 = header.find('\'', d + 8);
  const auto q2 = header.find('\'', q1 + 1);
  arr.dtype = header.substr(q1 + 1, q2 - q1 - 1);
  if (header.find("'fortran_order': True") != std::string::npos) {
    throw std::runtime_error("read_npy: Fortran order is not supported");
  }
  const auto s1 = header.find('(', header.find("'shape':"));
  const auto s2 = header.find(')', s1);
  std::string dims = header.substr(s1 + 1, s2 - s1 - 1);
  std::size_t count = 1;
  for (std::size_t p = 0; p < dims.size();) {
    const auto comma = dims.find(',', p);
    const std::string tok =
        dims.substr(p, comma == std::string::npos ? std::string::npos : comma - p);
    if (tok.find_first_of("0123456789") != std::string::npos) {
      arr.shape.push_back(static_cast<std::size_t>(std::stoull(tok)));
      count *= arr.shape.back();
    }
    if (comma == std::string::npos) {
      break;
    }
    p = comma + 1;
  }

  auto read_all = [&](auto tag, std::size_t n) {
    using T = decltype(tag);
    std::vector<T> buf(n);
    f.read(reinterpret_cast<char*>(buf.data()), static_cast<std::streamsize>(n * sizeof(T)));
    if (!f) {
      throw std::runtime_error("read_npy: truncated data in " + path.string());
    }
    return buf;
  };
  if (arr.dtype == "<f4") {
    const auto buf = read_all(float{}, count);
    arr.values.assign(buf.begin(), buf.end());
  } else if (arr.dtype == "<f8") {
    arr.values = read_all(double{}, count);
  } else if (arr.dtype == "|u1") {
    const auto buf = read_all(std::uint8_t{}, count);
    arr.values.assign(buf.begin(), buf.end());
  } else if (arr.dtype == "<c8") {
    const auto buf = read_all(std::complex<float>{}, count);
    arr.complex_data.assign(buf.begin(), buf.end());
    arr.values.resize(count);
    for (std::size_t i = 0; i < count; ++i) {
      arr.values[i] = std::abs(arr.complex_data[i]);
    }
  } else {
    throw std::runtime_error("read_npy: unsupported dtype " + arr.dtype);
  }
  return arr;
}

Image<std::uint8_t> to_db_u8(const Image<float>& linear, double min_db, double max_db) {
  Image<std::uint8_t> out(linear.width, linear.height);
  const double span = std::max(max_db - min_db, 1e-9);
  for (std::size_t i = 0; i < linear.size(); ++i) {
    const double v = linear.data[i];
    if (!(v > 0.0)) {
      continue;
    }
    const double t = (10.0 * std::log10(v) - min_db) / span;
    out.data[i] = static_cast<std::uint8_t>(std::lround(255.0 * std::clamp(t, 0.0, 1.0)));
  }
  return out;
}

Image<std::uint8_t> to_db_u8_auto(const Image<float>& linear, double dynamic_range_db) {
  std::vector<float> positive;
  for (const float v : linear.data) {
    if (v > 0.0f) {
      positive.push_back(v);
    }
  }
  if (positive.empty()) {
    return Image<std::uint8_t>(linear.width, linear.height);
  }
  auto percentile_db = [&](double q) {
    const auto k = static_cast<std::size_t>(q * static_cast<double>(positive.size() - 1));
    std::nth_element(positive.begin(), positive.begin() + static_cast<std::ptrdiff_t>(k),
                     positive.end());
    return 10.0 * std::log10(static_cast<double>(positive[k]));
  };
  // Saturate point-like scatterers (e.g. dihedrals 40+ dB above the
  // background) rather than letting them push the terrain to black.
  const double max_db = std::min(percentile_db(0.999), percentile_db(0.5) + 20.0);
  return to_db_u8(linear, max_db - dynamic_range_db, max_db);
}

Image<std::uint8_t> to_u8(const Image<float>& values, double lo, double hi) {
  Image<std::uint8_t> out(values.width, values.height);
  const double span = std::max(hi - lo, 1e-12);
  for (std::size_t i = 0; i < values.size(); ++i) {
    const double t = (static_cast<double>(values.data[i]) - lo) / span;
    out.data[i] = static_cast<std::uint8_t>(std::lround(255.0 * std::clamp(t, 0.0, 1.0)));
  }
  return out;
}

}  // namespace bsar
