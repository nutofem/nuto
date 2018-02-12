add_unit_test(UnstructuredGrid)
target_link_libraries(UnstructuredGrid Visualize)

add_unit_test(DataArray)

add_unit_test(VoronoiGeometries)
target_link_libraries(VoronoiGeometries Math)

add_unit_test(AverageHandler
    mechanics/elements/ElementShapeFunctions.cpp)
target_link_libraries(AverageHandler Visualize)

add_unit_test(VoronoiHandler
    mechanics/elements/ElementShapeFunctions.cpp)
target_link_libraries(VoronoiHandler Visualize)

add_unit_test(PointHandler
    mechanics/elements/ElementShapeFunctions.cpp
    mechanics/integrationtypes/IntegrationTypeBase.cpp
    mechanics/integrationtypes/IntegrationTypeTensorProduct.cpp)

target_link_libraries(PointHandler Visualize Math)

add_unit_test(Visualizer
    mechanics/elements/ElementShapeFunctions.cpp)
target_link_libraries(Visualizer Visualize)


add_unit_test(MeshFunction)
target_link_libraries(MeshFunction Mechanics Visualize)
