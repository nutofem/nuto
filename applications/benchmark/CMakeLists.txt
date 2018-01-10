unset(all_benchmarks CACHE)

add_benchmark(BuildGradient)

add_benchmark(ConstraintsBenchmark)

add_benchmark(DofMatrixExport)

add_benchmark(MeshFemBenchmark)

add_benchmark(NewtonRaphsonBenchmark)

add_benchmark(ShapeFunctionMemoization)

add_benchmark(VisualizeBenchmark)

add_custom_target(benchmarks DEPENDS ${all_benchmarks})
