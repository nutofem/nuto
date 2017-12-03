unset(all_benchmarks CACHE)

foreach(test
    BuildGradient
    LinearElasticity
    LinearElasticityExplicit
    LinearElasticityOPENMP
    LinearElasticityVisualize
    AssembleSparseMatrix
    ShapeFunctionMemoization
    Solvers
    MeshFemBenchmark
    )

    add_executable(${test} EXCLUDE_FROM_ALL ${test}.cpp)
    target_link_libraries(${test} Mechanics Math Base Visualize)
    append_to_benchmarks(${test})
endforeach()

if(PARDISO_FOUND)
    target_link_libraries(Solvers ${PARDISO_LIBRARIES} ${BLAS_LIBRARIES})
endif()

add_executable(DenseSpherePacking EXCLUDE_FROM_ALL DenseSpherePacking.cpp)
target_link_libraries(DenseSpherePacking Math Base GeometryConcrete)
add_test(benchmark::DenseSpherePacking
    ${CMAKE_CURRENT_BUILD_DIR}/DenseSpherePacking)
append_to_benchmarks(DenseSpherePacking)

add_executable(NewtonRaphsonBenchmark
    EXCLUDE_FROM_ALL NewtonRaphsonBenchmark.cpp)
target_include_directories(NewtonRaphsonBenchmark
    PRIVATE ${CMAKE_SOURCE_DIR}/src)
append_to_benchmarks(NewtonRaphsonBenchmark)

add_executable(BlockVectors
    EXCLUDE_FROM_ALL BlockVectors.cpp)
target_include_directories(BlockVectors
    PRIVATE ${CMAKE_SOURCE_DIR}/src)
append_to_benchmarks(BlockVectors)


add_executable(ConstraintsBenchmark
    EXCLUDE_FROM_ALL ConstraintsBenchmark.cpp
    ${CMAKE_SOURCE_DIR}/src/mechanics/constraintsPde/Constraints.cpp)
target_include_directories(ConstraintsBenchmark
    PRIVATE ${CMAKE_SOURCE_DIR}/src)
append_to_benchmarks(ConstraintsBenchmark)

add_executable(ParallelAssembly EXCLUDE_FROM_ALL ParallelAssembly.cpp)
target_link_libraries(ParallelAssembly
    Mechanics Visualize Boost::unit_test_framework)
target_include_directories(ParallelAssembly
    PRIVATE ${CMAKE_SOURCE_DIR}/test/tools)


add_benchmark(SolversGB)

add_custom_target(benchmarks DEPENDS ${all_benchmarks})
