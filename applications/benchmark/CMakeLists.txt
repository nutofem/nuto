add_custom_target(benchmark COMMAND ctest -V -R benchmark)


if(ENABLE_MECHANICS)
foreach(test
    BuildGradient
    LinearElasticity
    LinearElasticityVisualize
    AssembleSparseMatrix
    Solvers
    )

    add_executable(${test} ${test}.cpp)
    target_link_libraries(${test} NuToMechanics NuToMath NuToBase -lprofiler)
    if(MUMPS_FOUND)
        target_link_libraries(${test} -lprofiler ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
    endif()
    if(ENABLE_VISUALIZE)
        target_link_libraries(${test} NuToVisualize)
    endif()
    add_test(benchmark::${test} ${CMAKE_CURRENT_BINARY_DIR}/${test})
endforeach()

if(PARDISO_FOUND)
    target_link_libraries(Solvers ${PARDISO_LIBRARIES} ${OpenBLAS_LIBRARIES})
endif()


add_executable(LinearElasticityOPENMP LinearElasticityOPENMP.cpp)
target_link_libraries(LinearElasticityOPENMP NuToMechanics NuToMath NuToBase)
if(ENABLE_VISUALIZE)
    target_link_libraries(LinearElasticityOPENMP NuToVisualize)
endif()
endif()

if(ENABLE_GEOMETRYCONCRETE)
    add_executable(DenseSpherePacking DenseSpherePacking.cpp)
    target_link_libraries(DenseSpherePacking NuToMath NuToBase NuToGeometryConcrete ${Boost_LIBRARIES})
    add_test(benchmark::DenseSpherePacking ${CMAKE_CURRENT_BUILD_DIR}/DenseSpherePacking)
endif()
