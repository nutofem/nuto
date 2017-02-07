add_custom_target(benchmark COMMAND ctest -V -R benchmark)

# generate tests
foreach(test
    BuildGradient
    LinearElasticity
    AssembleSparseMatrix
    )

    add_executable(${test} ${test}.cpp)
    target_link_libraries(${test} NuToMechanics NuToMath NuToBase)
    if(MUMPS_FOUND)
        target_link_libraries(${test} ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
    endif()
    if(ENABLE_VISUALIZE)
        target_link_libraries(${test} NuToVisualize)
    endif()
    add_test(benchmark::${test} ${CMAKE_CURRENT_BINARY_DIR}/${test})
endforeach()

if(ENABLE_GEOMETRYCONCRETE)
    add_executable(DenseSpherePacking DenseSpherePacking.cpp)
    target_link_libraries(DenseSpherePacking NuToMath NuToBase NuToGeometryConcrete ${Boost_LIBRARIES})
    add_test(benchmark::DenseSpherePacking ${CMAKE_CURRENT_BUILD_DIR}/DenseSpherePacking)
endif()