#include "Benchmark.h"
#include "mechanics/mesh/MeshFemDofConvert.h"
#include "mechanics/mesh/UnitMeshFem.h"
#include "mechanics/interpolation/InterpolationTriangleLinear.h"

BENCHMARK(Mesh, CreateTriangles10x10, runner)
{
    while (runner.KeepRunningTime(1))
        NuTo::MeshFem mesh = NuTo::UnitMeshFem::CreateTriangles(10, 10);
}
BENCHMARK(Mesh, CreateTriangles100x100, runner)
{
    while (runner.KeepRunningTime(1))
        NuTo::MeshFem mesh = NuTo::UnitMeshFem::CreateTriangles(100, 100);
}
BENCHMARK(Mesh, CreateTriangles1000x1000, runner)
{
    while (runner.KeepRunningTime(1))
        NuTo::MeshFem mesh = NuTo::UnitMeshFem::CreateTriangles(1000, 1000);
}

BENCHMARK(Mesh, Convert10x10, runner)
{
    NuTo::InterpolationTriangleLinear interpolation(2);
    NuTo::MeshFem mesh = NuTo::UnitMeshFem::CreateTriangles(10, 10);
    while (runner.KeepRunningTime(1))
    {
        NuTo::MeshFem m = mesh;
        NuTo::AddDofInterpolation(&m, NuTo::DofType("dof", 2), interpolation);
    }
}
BENCHMARK(Mesh, Convert100x100, runner)
{
    NuTo::InterpolationTriangleLinear interpolation(2);
    NuTo::MeshFem mesh = NuTo::UnitMeshFem::CreateTriangles(100, 100);
    while (runner.KeepRunningTime(1))
    {
        NuTo::MeshFem m = mesh;
        NuTo::AddDofInterpolation(&m, NuTo::DofType("dof", 2), interpolation);
    }
}
BENCHMARK(Mesh, Convert1000x1000, runner)
{
    NuTo::InterpolationTriangleLinear interpolation(2);
    NuTo::MeshFem mesh = NuTo::UnitMeshFem::CreateTriangles(1000, 1000);
    while (runner.KeepRunningTime(1))
        NuTo::AddDofInterpolation(&mesh, NuTo::DofType("dof", 2), interpolation);
}


auto f = [](Eigen::VectorXd oldCoords) -> Eigen::VectorXd {
    return Eigen::Vector2d(oldCoords[0] * 12 + 1, oldCoords[1] * 4 - 6174);
};

BENCHMARK(Mesh, Transform10x10, runner)
{
    NuTo::MeshFem mesh = NuTo::UnitMeshFem::CreateTriangles(10, 10);
    while (runner.KeepRunningTime(1))
        NuTo::MeshFem newMesh = NuTo::UnitMeshFem::Transform(std::move(mesh), f);
}
BENCHMARK(Mesh, Transform100x100, runner)
{
    NuTo::MeshFem mesh = NuTo::UnitMeshFem::CreateTriangles(100, 100);
    while (runner.KeepRunningTime(1))
        NuTo::MeshFem newMesh = NuTo::UnitMeshFem::Transform(std::move(mesh), f);
}

BENCHMARK(Mesh, Transform1000x1000, runner)
{
    NuTo::MeshFem mesh = NuTo::UnitMeshFem::CreateTriangles(1000, 1000);
    while (runner.KeepRunningTime(1))
        NuTo::MeshFem newMesh = NuTo::UnitMeshFem::Transform(std::move(mesh), f);
}
