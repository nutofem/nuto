add_subdirectory(paraviewScripts)
add_subdirectory(meshes)

################################################################################
##  TESTS
################################################################################

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    foreach(test Feti FetiGradientDamage)
        add_executable(${test} EXCLUDE_FROM_ALL ${test}.cpp)
        target_link_libraries(${test} Mechanics Math Base FetiLib)
        target_link_libraries(${test} Boost::mpi ${LAPACK_LIBRARIES})
        target_link_libraries(${test} ${MPI_LIBRARIES})
        target_link_libraries(${test} Visualize)
        if(MUMPS_FOUND)
            target_link_libraries(${test} Mumps::Mumps)
            target_link_libraries(${test} ${CMAKE_THREAD_LIBS_INIT})
        endif()
        set(MPIEXEC_NUMPROC_FLAG -np 2)
        set(MPIEXEC mpirun)
        set(MPICMD ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${test})
        add_test(integration::mechanics::${test}.cpp ${MPICMD})
        append_to_tests(${test})
    endforeach()

################################################################################

    foreach(test FetiMultiplicityScaling FetiPreconditioner FetiKScaling)
        add_executable(${test} EXCLUDE_FROM_ALL ${test}.cpp)
        target_link_libraries(${test} Mechanics Math Base FetiLib)
        target_link_libraries(${test} Boost::mpi ${LAPACK_LIBRARIES})
        target_link_libraries(${test} ${MPI_LIBRARIES})
        target_link_libraries(${test} Visualize)
        set(MPIEXEC_NUMPROC_FLAG -np 4)
        set(MPIEXEC mpirun)
        set(MPICMD ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${test})
        add_test(integration::mechanics::${test}.cpp ${MPICMD})
        append_to_tests(${test})
    endforeach()

################################################################################

    foreach(test FetiRixenFarhatHeterogeneousExample)
        add_executable(${test} EXCLUDE_FROM_ALL ${test}.cpp)
        target_link_libraries(${test} Mechanics Math Base FetiLib)
        target_link_libraries(${test} Boost::mpi ${LAPACK_LIBRARIES})
        target_link_libraries(${test} ${MPI_LIBRARIES})
        target_link_libraries(${test} Visualize)
        set(MPIEXEC_NUMPROC_FLAG -np 16)
        set(MPIEXEC mpirun)
        set(MPICMD ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${test})
        add_test(integration::mechanics::${test}.cpp ${MPICMD})
        append_to_tests(${test})
    endforeach()
endif()
