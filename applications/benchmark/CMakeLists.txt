add_custom_target(benchmark COMMAND ctest -V -R benchmark)


if(ENABLE_MECHANICS)
foreach(test
    BuildGradient
    LinearElasticity
    LinearElasticityExplicit
    LinearElasticityVisualize
    AssembleSparseMatrix
    ShapeFunctionMemoization
    Solvers
    )

    add_executable(${test} ${test}.cpp)
    target_link_libraries(${test} Mechanics Math Base)
    if(MUMPS_FOUND)
        target_link_libraries(${test} ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
    endif()
    if(ENABLE_VISUALIZE)
        target_link_libraries(${test} Visualize)
    endif()
    add_test(benchmark::${test} ${CMAKE_CURRENT_BINARY_DIR}/${test})
endforeach()

if(PARDISO_FOUND)
    target_link_libraries(Solvers ${PARDISO_LIBRARIES} ${OpenBLAS_LIBRARIES})
endif()


add_executable(LinearElasticityOPENMP LinearElasticityOPENMP.cpp)
target_link_libraries(LinearElasticityOPENMP Mechanics Math Base)
if(ENABLE_VISUALIZE)
    target_link_libraries(LinearElasticityOPENMP Visualize)
endif()
endif()

if(ENABLE_GEOMETRYCONCRETE)
    add_executable(DenseSpherePacking DenseSpherePacking.cpp)
    target_link_libraries(DenseSpherePacking Math Base GeometryConcrete ${Boost_LIBRARIES})
    add_test(benchmark::DenseSpherePacking ${CMAKE_CURRENT_BUILD_DIR}/DenseSpherePacking)
endif()

add_executable(NewtonRaphsonBenchmark NewtonRaphsonBenchmark.cpp)
target_include_directories(NewtonRaphsonBenchmark PUBLIC ${CMAKE_SOURCE_DIR}/src)

