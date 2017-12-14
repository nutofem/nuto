#include <benchmark/benchmark.h>

#include "mechanics/mesh/MeshFemDofConvert.h"
#include "mechanics/mesh/UnitMeshFem.h"
#include "mechanics/interpolation/InterpolationTriangleLinear.h"

/*
 * Calculates the big O complexity of various, potentially expensive MeshFem methods.
 */

//! @brief Measures to create a 2D triangle mesh
static void Create(benchmark::State& state)
{
    const int n = state.range(0);
    for (auto _ : state)
        NuTo::MeshFem mesh = NuTo::UnitMeshFem::CreateTriangles(n, n);
    state.SetComplexityN(n * n);
}

//! @brief Measures to create a 2D triangle mesh and add another 'layer' of dof elements
static void Convert(benchmark::State& state)
{
    const int n = state.range(0);
    NuTo::InterpolationTriangleLinear interpolation(2);
    for (auto _ : state)
    {
        NuTo::MeshFem mesh = NuTo::UnitMeshFem::CreateTriangles(n, n);
        NuTo::AddDofInterpolation(&mesh, NuTo::DofType("dof", 2), interpolation);
    }
    state.SetComplexityN(n * n);
}

//! @brief Measures to transform the coordinates of a 2D triangle mesh
static void Transform(benchmark::State& state)
{
    const int n = state.range(0);

    auto f = [](Eigen::VectorXd oldCoords) -> Eigen::VectorXd {
        return Eigen::Vector2d(oldCoords[0] * 12 + 1, oldCoords[1] * 4 - 6174);
    };

    NuTo::MeshFem mesh = NuTo::UnitMeshFem::CreateTriangles(n, n);
    for (auto _ : state)
    {
        NuTo::MeshFem newMesh = NuTo::UnitMeshFem::Transform(std::move(mesh), f);
    }
    state.SetComplexityN(n * n);
}
BENCHMARK(Create)->RangeMultiplier(2)->Range(16, 1024)->Complexity();
BENCHMARK(Convert)->RangeMultiplier(2)->Range(16, 1024)->Complexity();
BENCHMARK(Transform)->RangeMultiplier(2)->Range(16, 1024)->Complexity();
BENCHMARK_MAIN()
