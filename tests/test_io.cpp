#include <catch2/catch_test_macros.hpp>
#include <filesystem>
#include <fstream>

#include "backscatter/image/io.hpp"

using namespace bsar;

namespace {
std::filesystem::path tmp(const char* name) {
  return std::filesystem::temp_directory_path() / name;
}
}  // namespace

TEST_CASE("NPY round trip", "[io]") {
  Image<float> f(5, 3);
  for (std::size_t i = 0; i < f.size(); ++i) {
    f.data[i] = static_cast<float>(i) * 0.5f;
  }
  write_npy(tmp("bsar_f.npy"), f);
  const NpyArray a = read_npy(tmp("bsar_f.npy"));
  CHECK(a.dtype == "<f4");
  CHECK(a.shape == std::vector<std::size_t>{3, 5});
  CHECK(a.values[7] == 3.5);

  Image<std::complex<float>> c(2, 2);
  c.data[3] = {3.0f, 4.0f};
  write_npy(tmp("bsar_c.npy"), c);
  const NpyArray b = read_npy(tmp("bsar_c.npy"));
  CHECK(b.dtype == "<c8");
  CHECK(b.values[3] == 5.0);
  CHECK(b.complex_data[3] == std::complex<double>(3.0, 4.0));

  // Header is padded to a multiple of 64 bytes.
  std::ifstream in(tmp("bsar_f.npy"), std::ios::binary);
  in.seekg(8);
  unsigned char len[2];
  in.read(reinterpret_cast<char*>(len), 2);
  CHECK((10 + len[0] + 256 * len[1]) % 64 == 0);
}

TEST_CASE("PNG writer produces a valid signature and size", "[io]") {
  Image<std::uint8_t> g(300, 300);
  for (std::size_t i = 0; i < g.size(); ++i) {
    g.data[i] = static_cast<std::uint8_t>(i % 251);
  }
  write_png(tmp("bsar_g.png"), g);  // > 65535 bytes: several stored blocks
  std::ifstream in(tmp("bsar_g.png"), std::ios::binary);
  char sig[8];
  in.read(sig, 8);
  CHECK(std::string(sig + 1, 3) == "PNG");
  CHECK(std::filesystem::file_size(tmp("bsar_g.png")) > 300u * 301u);
  CHECK_THROWS(write_png(tmp("bsar_bad.png"), 2, 2, 4, std::vector<std::uint8_t>(16)));
}

TEST_CASE("dB scaling", "[io]") {
  Image<float> f(3, 1);
  f.data = {0.0f, 0.01f, 1.0f};
  const auto u = to_db_u8(f, -20.0, 0.0);
  CHECK(u.data[0] == 0);
  CHECK(u.data[1] == 0);
  CHECK(u.data[2] == 255);
}
