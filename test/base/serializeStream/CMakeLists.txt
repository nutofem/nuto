add_executable(SerializeStream SerializeStream.cpp)
target_link_libraries(SerializeStream NuToMath NuToBase)
add_test(c::unit::SerializeStream ${CMAKE_CURRENT_BUILD_DIR}/SerializeStream)
