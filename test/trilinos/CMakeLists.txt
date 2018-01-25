#SET(Trilinos_PREFIX /home/briedel/install/Trilinos)
#SET(CMAKE_PREFIX_PATH ${Trilinos_PREFIX} ${CMAKE_PREFIX_PATH})
#FIND_PACKAGE(Trilinos REQUIRED)
#INCLUDE_DIRECTORIES(${Trilinos_INCLUDE_DIRS})
#INCLUDE_DIRECTORIES(${Trilinos_TPL_INCLUDE_DIRS})

SET(Trilinos_INCLUDE_DIRS /usr/include/trilinos)
SET(Trilinos_TPL_INCLUDE_DIRS /usr/include)

SET(MyTrilinos_LIBRARIES
    trilinos_epetra
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

INCLUDE_DIRECTORIES(${Trilinos_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${Trilinos_TPL_INCLUDE_DIRS})

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/test/tools)

#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/applications/custom)

#add_unit_test(TrilinosTests)
#target_link_libraries(TrilinosTests ${MPI_LIBRARIES} ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})



function(add_mpi_test name test_prefix num_proc bin_dir)
# Make sure to link MPI here too:
    target_link_libraries(${name} ${MPI_LIBRARIES})
    set(MPIEXEC_NUMPROC_FLAG -np ${num_proc})
    set(MPIEXEC mpirun)
    set(TEST_NAME ${test_prefix}::${name}.cpp)
    add_test(NAME ${TEST_NAME} COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${bin_dir}/${name})
endfunction(add_mpi_test)

#add_mpi_test(TrilinosTests 2)

foreach(test
        Conversion
        DirectSolver
        Epetra_Export
        Epetra_GraphCreation
        IterativeSolver)

    add_executable(${test} ${test}.cpp)
    #target_link_libraries(${test} Mechanics Math Base ConversionTools ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})
    target_link_libraries(${test} Mechanics Math Base ConversionTools Visualize ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${MyTrilinos_LIBRARIES})
    if(MUMPS_FOUND)
        target_link_libraries(${test} ${MUMPS_LIBRARIES})
    endif()

    set(NUM_PROC 2)
    add_mpi_test(${test} unit::trilinos ${NUM_PROC} ${CMAKE_CURRENT_BINARY_DIR})

endforeach()

