add_unit_test(UnstructuredGrid)
target_link_libraries(UnstructuredGrid NuTo)

add_unit_test(DataArray)

add_unit_test(VoronoiGeometries)
target_link_libraries(VoronoiGeometries NuTo)

add_unit_test(AverageHandler)
target_link_libraries(AverageHandler NuTo)

add_unit_test(VoronoiHandler)
target_link_libraries(VoronoiHandler NuTo)

add_unit_test(PointHandler)
target_link_libraries(PointHandler NuTo)

add_unit_test(Visualizer)
target_link_libraries(Visualizer NuTo)

add_unit_test(MeshFunction)
target_link_libraries(MeshFunction NuTo)

add_unit_test(PostProcess)
target_link_libraries(PostProcess NuTo Boost::filesystem)

