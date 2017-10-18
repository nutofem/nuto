pybind11_add_module(nuto HygroHelpers.cpp ${CMAKE_SOURCE_DIR}/src/mechanics/PDEs/HygroHelpers.cpp)
target_include_directories(nuto PRIVATE ${CMAKE_SOURCE_DIR}/src)
