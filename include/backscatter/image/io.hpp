#pragma once

#include <complex>
#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

#include "backscatter/image/image.hpp"

namespace bsar {

/// 8-bit grayscale (channels = 1) or RGB (channels = 3) PNG. Uses stored
/// (uncompressed) deflate blocks, so no zlib dependency.
void write_png(const std::filesystem::path& path, std::size_t width, std::size_t height,
               int channels, const std::vector<std::uint8_t>& pixels);
void write_png(const std::filesystem::path& path, const Image<std::uint8_t>& gray);

/// NumPy .npy files (format 1.0, little-endian, C order, shape (height, width)).
void write_npy(const std::filesystem::path& path, const Image<float>& image);
void write_npy(const std::filesystem::path& path, const Image<std::complex<float>>& image);
void write_npy(const std::filesystem::path& path, const Image<std::uint8_t>& image);

/// A loaded .npy array of dtype <f4, <f8, <c8 or |u1, converted to doubles
/// (complex arrays are converted to magnitude unless `complex_data` is read).
struct NpyArray {
  std::vector<std::size_t> shape;
  std::string dtype;
  std::vector<double> values;                      // real data or magnitudes
  std::vector<std::complex<double>> complex_data;  // only for <c8
};
NpyArray read_npy(const std::filesystem::path& path);

/// Map values to 8 bits on a decibel scale, clamping to [min_db, max_db].
/// Non-positive inputs map to 0.
Image<std::uint8_t> to_db_u8(const Image<float>& linear, double min_db, double max_db);
/// Same, choosing the range from the data: the top is the 99.9th percentile of
/// positive values, capped at the median + 20 dB.
Image<std::uint8_t> to_db_u8_auto(const Image<float>& linear, double dynamic_range_db = 40.0);
/// Linear map of [lo, hi] to [0, 255].
Image<std::uint8_t> to_u8(const Image<float>& values, double lo, double hi);

}  // namespace bsar
