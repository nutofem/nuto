create_symlink("MeshCompanionGmsh.msh")

add_executable(MeshCompanion MeshCompanion.cpp)
target_link_libraries(MeshCompanion NuToMechanics NuToMath NuToBase ${Boost_LIBRARIES})
if(ENABLE_VISUALIZE)
    target_link_libraries(MeshCompanion NuToVisualize)
endif()
add_test(almost::mechanics::mesh::MeshCompanion ${CMAKE_CURRENT_BINARY_DIR}/MeshCompanion)

add_unit_test(Mesh
    mechanics/elements/ElementShapeFunctions)
