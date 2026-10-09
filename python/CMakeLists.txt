find_package(Python 3.9 COMPONENTS Interpreter Development.Module REQUIRED)

# Locate nanobind from the active Python environment (pip install nanobind).
execute_process(
  COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
  OUTPUT_STRIP_TRAILING_WHITESPACE
  OUTPUT_VARIABLE nanobind_ROOT)
find_package(nanobind CONFIG REQUIRED)

nanobind_add_module(backscatter_py NB_STATIC backscatter_module.cpp)
set_target_properties(backscatter_py PROPERTIES OUTPUT_NAME backscatter)
target_link_libraries(backscatter_py PRIVATE backscatter::backscatter)

# The static core library is linked into a shared module.
set_target_properties(backscatter PROPERTIES POSITION_INDEPENDENT_CODE ON)
