add_executable(bsar_bench bsar_bench.cpp)
target_link_libraries(bsar_bench PRIVATE backscatter::backscatter)
bsar_set_warnings(bsar_bench)

# Intel Embree is a benchmark baseline only, never a dependency of the core.
if(BSAR_WITH_EMBREE)
  find_package(embree 4 CONFIG REQUIRED)
  target_link_libraries(bsar_bench PRIVATE embree)
  target_compile_definitions(bsar_bench PRIVATE BSAR_HAVE_EMBREE)
endif()
