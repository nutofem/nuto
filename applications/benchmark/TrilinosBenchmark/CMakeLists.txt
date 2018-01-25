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

if(ENABLE_MECHANICS)
foreach(test
    IterativeSolvers
    DirectSolvers
    Conversion_NuTo_Trilinos
    )

    add_executable(${test} ${test}.cpp)
    #target_link_libraries(${test} Mechanics Math Base ConversionTools ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})
    target_link_libraries(${test} Mechanics Math Base ConversionTools ${MyTrilinos_LIBRARIES})
    if(ENABLE_VISUALIZE)
        target_link_libraries(${test} Visualize)
    endif()
    add_test(NAME benchmark::Trilinos::${test}.cpp COMMAND mpirun -np 4 ${CMAKE_CURRENT_BINARY_DIR}/${test})
endforeach()

endif()
