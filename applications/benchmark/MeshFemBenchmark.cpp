#include <benchmark/benchmark.h>
#include <chrono>

#include "mechanics/mesh/MeshFemDofConvert.h"
#include "mechanics/mesh/UnitMeshFem.h"
#include "mechanics/interpolation/InterpolationTriangleLinear.h"

/*
 * Calculates the big O complexity of various, potentially expensive MeshFem methods.
 *
 * It may be a good example on how to use manual timing. Not all work done inside the
 * benchmarks loops is supposed to contribute to the big O complexity. These parts
 * are excluded from the measurement.
 */

//! @brief Measures time to create a 2D triangle mesh
static void Create(benchmark::State& state)
{
    const int n = state.range(0);
    for (auto _ : state)
        NuTo::MeshFem mesh = NuTo::UnitMeshFem::CreateTriangles(n, n);
    state.SetComplexityN(n * n);
}

//! @brief Measures time add another 'layer' of dof elements
static void Convert(benchmark::State& state)
{
    const int n = state.range(0);
    NuTo::InterpolationTriangleLinear interpolation;
    for (auto _ : state)
    {
        NuTo::MeshFem mesh = NuTo::UnitMeshFem::CreateTriangles(n, n);
        auto start = std::chrono::high_resolution_clock::now();
        NuTo::AddDofInterpolation(&mesh, NuTo::DofType("dof", 2), interpolation);
        auto end = std::chrono::high_resolution_clock::now();

        auto elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
        state.SetIterationTime(elapsed_seconds.count());
    }
    state.SetComplexityN(n * n);
}

//! @brief Measures time to transform the coordinates of a 2D triangle mesh
static void Transform(benchmark::State& state)
{
    const int n = state.range(0);

    auto f = [](Eigen::VectorXd oldCoords) -> Eigen::VectorXd {
        return Eigen::Vector2d(oldCoords[0] * 12 + 1, oldCoords[1] * 4 - 6174);
    };

    for (auto _ : state)
    {
        NuTo::MeshFem mesh = NuTo::UnitMeshFem::CreateTriangles(n, n);
        auto start = std::chrono::high_resolution_clock::now();
        NuTo::MeshFem newMesh = NuTo::UnitMeshFem::Transform(std::move(mesh), f);
        auto end = std::chrono::high_resolution_clock::now();

        auto elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
        state.SetIterationTime(elapsed_seconds.count());
    }
    state.SetComplexityN(n * n);
}
BENCHMARK(Create)->RangeMultiplier(2)->Range(16, 1024)->Complexity();
BENCHMARK(Convert)->RangeMultiplier(2)->Range(16, 1024)->Complexity();
BENCHMARK(Transform)->RangeMultiplier(2)->Range(16, 1024)->Complexity();
BENCHMARK_MAIN();
