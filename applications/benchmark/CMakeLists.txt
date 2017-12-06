unset(all_benchmarks CACHE)

add_benchmark(ShapeFunctionMemoization)

add_benchmark(MeshFemBenchmark)

add_benchmark(DenseSpherePacking)
target_link_libraries(DenseSpherePacking GeometryConcrete)

add_benchmark(NewtonRaphsonBenchmark)

add_benchmark(ConstraintsBenchmark)

#add_benchmark(SolversGB)

add_custom_target(benchmarks DEPENDS ${all_benchmarks})
