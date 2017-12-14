pybind11_add_module(hygro_helpers HygroHelpers.cpp ${CMAKE_SOURCE_DIR}/src/mechanics/PDEs/HygroHelpers.cpp)
target_include_directories(hygro_helpers PRIVATE ${CMAKE_SOURCE_DIR}/src)
