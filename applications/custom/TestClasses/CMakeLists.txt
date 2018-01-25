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


SET(Trilinos_useInstalledLibs False)

if (${Trilinos_useInstalledLibs})
    SET(Trilinos_PREFIX /home/briedel/install/Trilinos)
    SET(CMAKE_PREFIX_PATH ${Trilinos_PREFIX} ${CMAKE_PREFIX_PATH})
    FIND_PACKAGE(Trilinos REQUIRED)
else()
    SET(Trilinos_PREFIX /usr/include/trilinos)
    SET(Trilinos_INCLUDE_DIRS /usr/include/trilinos)
    SET(Trilinos_TPL_INCLUDE_DIRS /usr/include)
    SET(Trilinos_LIBRARIES /usr/lib/*.a)

    SET(MyTrilinos_LIBRARIES
        trilinos_epetra
        trilinos_epetraext
        trilinos_tpetra
        trilinos_teuchoscomm
        trilinos_teuchoscore
        trilinos_teuchoskokkoscomm
        trilinos_teuchoskokkoscompat
        trilinos_teuchosnumerics
        trilinos_teuchosparameterlist
        trilinos_kokkoscore
        trilinos_amesos2
        trilinos_amesos
        trilinos_aztecoo
        trilinos_ifpack2
        trilinos_ifpack
        trilinos_belos
        trilinos_belosepetra
        trilinos_belostpetra
        )
endif()

INCLUDE_DIRECTORIES(${Trilinos_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${Trilinos_TPL_INCLUDE_DIRS})

foreach(myTestFiles
    nuto_trilinos
    1D_Test
    2D_Test
    MeshJson_Test)
#    LoadBalancing)

    target_link_libraries(${myTestFiles} Mechanics Math Visualize ${MPI_Libraries})
    #target_link_libraries(${myTestFiles} ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})
    target_link_libraries(${myTestFiles} ${MyTrilinos_LIBRARIES} ${Boost_LIBRARIES})

endforeach()

add_library(StructureMesh STATIC StructureMesh.cpp)
#target_link_libraries(StructureMesh Mechanics Math Visualize ${Trilinos_LIBRARIES})
target_link_libraries(StructureMesh Mechanics Math Visualize ${Boost_LIBRARIES})

add_library(ConversionTools STATIC ConversionTools.cpp)
#target_link_libraries(ConversionTools ${Trilinos_LIBRARIES})
add_library(PrintTools PrintTools.cpp)
add_library(TrilinosUtils TrilinosUtils.cpp)
target_link_libraries(1D_Test ConversionTools PrintTools)
target_link_libraries(2D_Test ConversionTools PrintTools TrilinosUtils)

target_link_libraries(MeshJson_Test ConversionTools PrintTools FetiLib)

add_library(ZoltanMesh ZoltanMesh.cpp)
target_link_libraries(ZoltanMesh Mechanics Math Visualize)
