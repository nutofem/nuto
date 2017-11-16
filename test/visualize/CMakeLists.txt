
add_unit_test(UnstructuredGrid)
target_link_libraries(UnstructuredGrid Visualize)

add_unit_test(DataArray)

add_unit_test(GroupAverage
    mechanics/elements/ElementShapeFunctions.cpp)
target_link_libraries(GroupAverage Visualize)

add_unit_test(GroupVoronoiQuad
    mechanics/elements/ElementShapeFunctions.cpp)
target_link_libraries(GroupVoronoiQuad Visualize)
