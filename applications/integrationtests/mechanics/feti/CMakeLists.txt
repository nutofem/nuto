
########################################################################################################################
##  MESH FILES
########################################################################################################################

set(MESH_DIR "${CMAKE_CURRENT_SOURCE_DIR}/meshes")
foreach(rank RANGE 0 3)
    foreach(mesh FetiMultiplicityScaling.mesh)
        execute_process(COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink" "${MESH_DIR}/${mesh}${rank}" "${CMAKE_CURRENT_BINARY_DIR}/${mesh}${rank}")
    endforeach()
endforeach()

foreach(rank RANGE 0 1)
    foreach(mesh Feti.mesh)
        execute_process(COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink" "${MESH_DIR}/${mesh}${rank}" "${CMAKE_CURRENT_BINARY_DIR}/${mesh}${rank}")
    endforeach()
endforeach()

foreach(rank RANGE 0 15)
    foreach(mesh FetiRixenFarhatHeterogeneousExample.mesh)
        execute_process(COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink" "${MESH_DIR}/${mesh}${rank}" "${CMAKE_CURRENT_BINARY_DIR}/${mesh}${rank}")
    endforeach()
endforeach()

execute_process(COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink" "${MESH_DIR}/Feti_compare.msh" "${CMAKE_CURRENT_BINARY_DIR}/Feti_compare.msh")

########################################################################################################################
##  TESTS
########################################################################################################################

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    foreach(test Feti FetiGradientDamage)
        add_executable(${test} ${test}.cpp)
        target_link_libraries(${test} NuToMechanics NuToMath NuToBase NuToFeti ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${ANN_LIBRARIES} ${MPI_LIBRARIES})
        if(ENABLE_VISUALIZE)
            target_link_libraries(${test} NuToVisualize)
        endif()
        if(MUMPS_FOUND)
            target_link_libraries(${test} ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
        endif()
        set(MPIEXEC_NUMPROC_FLAG -np 2)
        set(MPIEXEC mpirun)
        add_test(NAME integration::mechanics::${test}.cpp COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${CMAKE_CURRENT_BINARY_DIR}/${test})
    endforeach()
    foreach(test FetiMultiplicityScaling FetiPreconditioners)
        add_executable(${test} ${test}.cpp)
        target_link_libraries(${test} NuToMechanics NuToMath NuToBase NuToFeti ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${ANN_LIBRARIES} ${MPI_LIBRARIES})
        if(ENABLE_VISUALIZE)
            target_link_libraries(${test} NuToVisualize)
        endif()
        set(MPIEXEC_NUMPROC_FLAG -np 4)
        set(MPIEXEC mpirun)
        add_test(NAME integration::mechanics::${test}.cpp COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${CMAKE_CURRENT_BINARY_DIR}/${test})
    endforeach()

    foreach(test FetiRixenFarhatHeterogeneousExample)
        add_executable(${test} ${test}.cpp)
        target_link_libraries(${test} NuToMechanics NuToMath NuToBase NuToFeti ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${ANN_LIBRARIES} ${MPI_LIBRARIES})
        if(ENABLE_VISUALIZE)
            target_link_libraries(${test} NuToVisualize)
        endif()
        set(MPIEXEC_NUMPROC_FLAG -np 16)
        set(MPIEXEC mpirun)
        add_test(NAME integration::mechanics::${test}.cpp COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${CMAKE_CURRENT_BINARY_DIR}/${test})
    endforeach()

endif()
