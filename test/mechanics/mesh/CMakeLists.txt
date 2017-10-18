create_symlink("MeshCompanionGmsh.msh")
create_symlink("MeshCompanionGmshFine.msh")

add_executable(MeshCompanion MeshCompanion.cpp)
target_link_libraries(MeshCompanion Mechanics Math Base ${Boost_LIBRARIES})
if(ENABLE_VISUALIZE)
    target_link_libraries(MeshCompanion Visualize)
endif()
add_test(almost::mechanics::mesh::MeshCompanion ${CMAKE_CURRENT_BINARY_DIR}/MeshCompanion)

add_unit_test(MeshFem
    mechanics/elements/ElementShapeFunctions)
target_link_libraries(MeshFem ${ANN_LIBRARIES})
