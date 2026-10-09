add_library(
  backscatter
  math/fft.cpp
  geometry/mesh.cpp
  geometry/bvh.cpp
  geometry/wide_bvh.cpp
  scene/material.cpp
  scene/heightmap.cpp
  scene/buildings.cpp
  scene/scene.cpp
  sensor/geodesy.cpp
  sensor/orbit.cpp
  sensor/radar.cpp
  sensor/sar_geometry.cpp
  trace/geometric.cpp
  trace/hillshade.cpp
  trace/coherent.cpp
  image/io.cpp
  image/backprojection.cpp
  image/backprojection_simd.cpp
  util/cpu.cpp
  eval/metrics.cpp
  eval/theory_checks.cpp)
add_library(backscatter::backscatter ALIAS backscatter)

target_include_directories(backscatter PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>)
target_compile_features(backscatter PUBLIC cxx_std_20)
target_link_libraries(backscatter PUBLIC Threads::Threads)
bsar_set_warnings(backscatter)

# The watertight ray/triangle test relies on a shared edge's 2D edge function
# being computed as exact negatives by the two triangles that share it. FMA
# contraction of `a*b - c*d` breaks that antisymmetry and lets grazing rays leak
# through seams (Clang contracts by default; GCC does outside strict ISO mode).
# PUBLIC because the test is inlined into consumers through triangle.hpp.
target_compile_options(backscatter PUBLIC $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-ffp-contract=off>)

if(BSAR_WITH_GDAL)
  find_package(GDAL CONFIG)
  if(NOT GDAL_FOUND)
    find_package(GDAL REQUIRED)
  endif()
  target_link_libraries(backscatter PUBLIC GDAL::GDAL)
  target_compile_definitions(backscatter PUBLIC BSAR_HAVE_GDAL)
endif()
