find_package(Catch2 3 QUIET)
if(NOT Catch2_FOUND)
  include(FetchContent)
  FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG v3.7.1
    GIT_SHALLOW TRUE)
  FetchContent_MakeAvailable(Catch2)
  list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
endif()
include(Catch)

add_executable(
  bsar_tests
  test_math.cpp
  test_bvh.cpp
  test_sensor.cpp
  test_scene.cpp
  test_geometric.cpp
  test_coherent.cpp
  test_eval.cpp
  test_io.cpp)
target_link_libraries(bsar_tests PRIVATE backscatter::backscatter Catch2::Catch2WithMain)
bsar_set_warnings(bsar_tests)

catch_discover_tests(bsar_tests DISCOVERY_MODE PRE_TEST)
