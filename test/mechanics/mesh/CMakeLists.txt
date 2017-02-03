execute_process(COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink" "${CMAKE_CURRENT_SOURCE_DIR}/MeshCompanionGmsh.msh" "${CMAKE_CURRENT_BINARY_DIR}/MeshCompanionGmsh.msh")

add_executable(MeshCompanion MeshCompanion.cpp)
target_link_libraries(MeshCompanion NuToMechanics NuToMath NuToBase ${Boost_LIBRARIES})
if(ENABLE_VISUALIZE)
    target_link_libraries(MeshCompanion NuToVisualize)
endif()
add_test(unit::mechanics::mesh::MeshCompanion ${CMAKE_CURRENT_BINARY_DIR}/MeshCompanion)
