add_unit_test(UnstructuredGrid)
target_link_libraries(UnstructuredGrid Visualize)

add_unit_test(DataArray)

add_unit_test(AverageVisualizer
    mechanics/elements/ElementShapeFunctions.cpp)
target_link_libraries(AverageVisualizer Visualize)

add_unit_test(VoronoiVisualizer
    mechanics/elements/ElementShapeFunctions.cpp)
target_link_libraries(VoronoiVisualizer Visualize)
