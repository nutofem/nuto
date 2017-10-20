# copy meshes
foreach(rank
    _000001
    _000002
    )
    execute_process( COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink" "${CMAKE_SOURCE_DIR}/applications/integrationtests/mechanics/meshes/trilinos.msh${rank}" "${CMAKE_CURRENT_BINARY_DIR}/trilinos.msh${rank}")
endforeach()


set(MESH_DIR "${CMAKE_CURRENT_SOURCE_DIR}/meshes")
foreach(mesh jsonMesh_Test)
    execute_process(COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink" "${MESH_DIR}/${mesh}.mesh" "${CMAKE_CURRENT_BINARY_DIR}/meshes/${mesh}.mesh")
    foreach(rank RANGE 0 1)
        execute_process(COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink" "${MESH_DIR}/${mesh}_0${rank}.mesh" "${CMAKE_CURRENT_BINARY_DIR}/meshes/${mesh}_0${rank}.mesh")
    endforeach()
endforeach()

foreach(myTestFiles
    nuto_trilinos
    1D_Test
    2D_Test
    MeshJson_Test)
#    LoadBalancing)

    add_executable(${myTestFiles} ${myTestFiles}.cpp)
#    target_link_libraries(${myTestFiles} NuToMechanics NuToVisualize ${MPI_Libraries})

endforeach()


SET(Trilinos_PREFIX /home/briedel/install/Trilinos)
SET(CMAKE_PREFIX_PATH ${Trilinos_PREFIX} ${CMAKE_PREFIX_PATH})
FIND_PACKAGE(Trilinos REQUIRED)
INCLUDE_DIRECTORIES(${Trilinos_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${Trilinos_TPL_INCLUDE_DIRS})


foreach(myTestFiles
    nuto_trilinos
    1D_Test
    2D_Test
    MeshJson_Test)
#    LoadBalancing)

    target_link_libraries(${myTestFiles} Mechanics Math Visualize ${MPI_Libraries})
    target_link_libraries(${myTestFiles} ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})

endforeach()

#add_library(StructureMesh StructureMesh.cpp)
#target_link_libraries(StructureMesh Mechanics Math Visualize)

add_library(ConversionTools ConversionTools.cpp)
add_library(PrintTools PrintTools.cpp)
add_library(TrilinosUtils TrilinosUtils.cpp)
target_link_libraries(1D_Test ConversionTools PrintTools)
target_link_libraries(2D_Test ConversionTools PrintTools TrilinosUtils)

target_link_libraries(MeshJson_Test ConversionTools PrintTools FetiLib)
