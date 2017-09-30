add_subdirectory(paraviewScripts)
add_subdirectory(meshes)

################################################################################
##  TESTS
################################################################################

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    foreach(test Feti FetiGradientDamage)
        add_executable(${test} ${test}.cpp)
        target_link_libraries(${test} Mechanics Math Base FetiLib)
        target_link_libraries(${test} Boost::mpi ${LAPACK_LIBRARIES})
        target_link_libraries(${test} ${ANN_LIBRARIES} ${MPI_LIBRARIES})
        if(ENABLE_VISUALIZE)
            target_link_libraries(${test} Visualize)
        endif()
        if(MUMPS_FOUND)
            target_link_libraries(${test} ${MUMPS_LIBRARIES})
            target_link_libraries(${test} ${CMAKE_THREAD_LIBS_INIT})
        endif()
        set(MPIEXEC_NUMPROC_FLAG -np 2)
        set(MPIEXEC mpirun)
        set(MPICMD ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${test})
        add_test(NAME integration::mechanics::${test}.cpp COMMAND ${MPICMD})
    endforeach()

################################################################################

    foreach(test FetiMultiplicityScaling FetiPreconditioner FetiKScaling)
        add_executable(${test} ${test}.cpp)
        target_link_libraries(${test} Mechanics Math Base FetiLib)
        target_link_libraries(${test} Boost::mpi ${LAPACK_LIBRARIES})
        target_link_libraries(${test} ${ANN_LIBRARIES} ${MPI_LIBRARIES})
        if(ENABLE_VISUALIZE)
            target_link_libraries(${test} Visualize)
        endif()
        set(MPIEXEC_NUMPROC_FLAG -np 4)
        set(MPIEXEC mpirun)
        set(MPICMD ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${test})
        add_test(NAME integration::mechanics::${test}.cpp COMMAND ${MPICMD})
    endforeach()

################################################################################

    foreach(test FetiRixenFarhatHeterogeneousExample)
        add_executable(${test} ${test}.cpp)
        target_link_libraries(${test} Mechanics Math Base FetiLib)
        target_link_libraries(${test} Boost::mpi ${LAPACK_LIBRARIES})
        target_link_libraries(${test} ${ANN_LIBRARIES} ${MPI_LIBRARIES})
        if(ENABLE_VISUALIZE)
            target_link_libraries(${test} Visualize)
        endif()
        set(MPIEXEC_NUMPROC_FLAG -np 16)
        set(MPIEXEC mpirun)
        set(MPICMD ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${test})
        add_test(NAME integration::mechanics::${test}.cpp COMMAND ${MPICMD})
    endforeach()
endif()
