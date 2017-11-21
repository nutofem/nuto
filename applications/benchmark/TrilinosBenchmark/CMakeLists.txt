SET(Trilinos_PREFIX /home/briedel/install/Trilinos)
SET(CMAKE_PREFIX_PATH ${Trilinos_PREFIX} ${CMAKE_PREFIX_PATH})
FIND_PACKAGE(Trilinos REQUIRED)
INCLUDE_DIRECTORIES(${Trilinos_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${Trilinos_TPL_INCLUDE_DIRS})

if(ENABLE_MECHANICS)
foreach(test
    IterativeSolvers
    DirectSolvers
    Conversion_NuTo_Trilinos
    )

    add_executable(${test} ${test}.cpp)
    target_link_libraries(${test} Mechanics Math Base ConversionTools ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})
    if(ENABLE_VISUALIZE)
        target_link_libraries(${test} Visualize)
    endif()
    add_test(NAME benchmark::Trilinos::${test}.cpp COMMAND mpirun -np 4 ${CMAKE_CURRENT_BINARY_DIR}/${test})
endforeach()

endif()
