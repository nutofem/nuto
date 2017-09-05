add_executable(MixedVisualization MixedVisualization.cpp)
target_link_libraries(MixedVisualization Math Base Visualize)
add_test(integration::visualize::MixedVisualization.cpp
    ${CMAKE_CURRENT_BUILD_DIR}/MixedVisualization)
