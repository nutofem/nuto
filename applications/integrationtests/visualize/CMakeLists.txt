add_executable(MixedVisualization EXCLUDE_FROM_ALL MixedVisualization.cpp)
target_link_libraries(MixedVisualization Math Base Visualize
    Boost::unit_test_framework)
add_test(integration::visualize::MixedVisualization.cpp
    ${CMAKE_CURRENT_BUILD_DIR}/MixedVisualization)
append_to_tests(MixedVisualization)
