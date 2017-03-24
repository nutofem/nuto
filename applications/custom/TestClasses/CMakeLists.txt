# copy meshes
foreach(rank
    _000001
    _000002
    )
    execute_process( COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink" "${CMAKE_SOURCE_DIR}/applications/integrationtests/mechanics/meshes/trilinos.msh${rank}" "${CMAKE_CURRENT_BINARY_DIR}/trilinos.msh${rank}")
endforeach()


foreach(myTestFiles
    nuto_trilinos)
    add_executable(${myTestFiles} ${myTestFiles}.cpp)
    target_link_libraries(${myTestFiles} NuToMechanics NuToMath NuToVisualize ${MPI_Libraries})

endforeach()


SET(Trilinos_PREFIX /home/briedel/install/Trilinos)
SET(CMAKE_PREFIX_PATH ${Trilinos_PREFIX} ${CMAKE_PREFIX_PATH})
FIND_PACKAGE(Trilinos REQUIRED)
INCLUDE_DIRECTORIES(${Trilinos_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${Trilinos_TPL_INCLUDE_DIRS})
#LINK_DIRECTORIES(${Trilinos_LIBRARY_DIRS})
#LINK_DIRECTORIES(${Trilinos_TPL_LIBRARY_DIRS})


foreach(myTestFiles
    nuto_trilinos)
    target_link_libraries(${myTestFiles} ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})

endforeach()

