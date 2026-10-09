add_executable(bsar_bench bsar_bench.cpp)
target_link_libraries(bsar_bench PRIVATE backscatter::backscatter)
bsar_set_warnings(bsar_bench)
