cmake_minimum_required(VERSION 3.21)

project(
  backscatter
  VERSION 0.1.0
  DESCRIPTION "A physically based SAR ray tracer"
  LANGUAGES CXX)

option(BSAR_BUILD_TESTS "Build the unit tests" ${PROJECT_IS_TOP_LEVEL})
option(BSAR_BUILD_APPS "Build the command-line tools" ${PROJECT_IS_TOP_LEVEL})
option(BSAR_BUILD_BENCH "Build the benchmark harness" ${PROJECT_IS_TOP_LEVEL})
option(BSAR_BUILD_PYTHON "Build the nanobind Python module" OFF)
option(BSAR_WITH_GDAL "Enable GeoTIFF DEM loading through GDAL" OFF)
option(BSAR_ENABLE_SANITIZERS "Build everything with ASan + UBSan" OFF)
option(BSAR_WARNINGS_AS_ERRORS "Treat compiler warnings as errors" OFF)
option(BSAR_NATIVE "Compile for the host CPU (-march=native)" OFF)
option(BSAR_WITH_EMBREE "Add Intel Embree 4 as a baseline in bsar_bench" OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include(Warnings)
include(Sanitizers)

if(BSAR_NATIVE AND NOT MSVC)
  add_compile_options(-march=native)
endif()

find_package(Threads REQUIRED)

add_subdirectory(src)

if(BSAR_BUILD_APPS)
  add_subdirectory(apps)
endif()

if(BSAR_BUILD_BENCH)
  add_subdirectory(bench)
endif()

if(BSAR_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()

if(BSAR_BUILD_PYTHON)
  add_subdirectory(python)
endif()
