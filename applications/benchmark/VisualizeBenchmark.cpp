#include <benchmark/benchmark.h>
#include "mechanics/mesh/UnitMeshFem.h"
#include "mechanics/mesh/MeshFemDofConvert.h"
#include "mechanics/cell/Cell.h"
#include "mechanics/interpolation/InterpolationQuadLinear.h"
#include "mechanics/interpolation/InterpolationQuadSerendipity.h"
#include "mechanics/integrationtypes/IntegrationTypeTensorProduct.h"

#include "visualize/Visualizer.h"
#include "visualize/VoronoiHandler.h"
#include "visualize/VoronoiGeometries.h"

/*
 * Find potential bottlenecks in visualization
 */


//! contains a mesh with num x num quad elements, linear dof interpolation and 4 integration points per cell
struct TestCells
{
    TestCells(int num)
        : mesh(NuTo::UnitMeshFem::CreateQuads(num, num))
    {
        const auto& interpolation = mesh.CreateInterpolation(NuTo::InterpolationQuadLinear());
        NuTo::AddDofInterpolation(&mesh, dof, interpolation);

        for (size_t id = 0; id < mesh.Elements.Size(); ++id)
            cells.push_back(NuTo::Cell(mesh.Elements[id], integrationType, id));
        for (auto& cell : cells)
            cellInterfaces.Add(cell);
    }

    NuTo::DofType dof = NuTo::DofType("displ", 2);
    NuTo::MeshFem mesh;
    NuTo::IntegrationTypeTensorProduct<2> integrationType =
            NuTo::IntegrationTypeTensorProduct<2>(2, NuTo::eIntegrationMethod::GAUSS);

    std::vector<NuTo::Cell> cells;
    NuTo::Group<NuTo::CellInterface> cellInterfaces;
};

//! benchmarks the creation of the geometry (creating points and cells)
void VisualizeQuadGeometry(benchmark::State& state)
{
    using namespace NuTo::Visualize;
    TestCells s(state.range(0));
    for (auto _ : state)
        Visualizer<VoronoiHandler> visu(s.cellInterfaces, VoronoiGeometryQuad(2));

    state.SetComplexityN(state.range(0) * state.range(0));
}
BENCHMARK(VisualizeQuadGeometry)->RangeMultiplier(2)->Range(2, 512)->Complexity();


//! benchmarks the extraction and writing of point data
void VisualizeQuadPointData(benchmark::State& state)
{
    using namespace NuTo::Visualize;
    TestCells s(state.range(0));

    Visualizer<VoronoiHandler> visu(s.cellInterfaces, VoronoiGeometryQuad(2));

    for (auto _ : state)
        visu.DofValues(s.dof);

    state.SetComplexityN(state.range(0) * state.range(0));
}
BENCHMARK(VisualizeQuadPointData)->RangeMultiplier(2)->Range(2, 512)->Complexity();


//! benchmarks the extraction and writing of cell data
void VisualizeQuadCellData(benchmark::State& state)
{
    using namespace NuTo::Visualize;
    TestCells s(state.range(0));

    Visualizer<VoronoiHandler> visu(s.cellInterfaces, VoronoiGeometryQuad(2));
    auto f = [&](const NuTo::CellData&, const NuTo::CellIpData&) { return Eigen::Vector2d(1, 2); };

    for (auto _ : state)
        visu.CellData(f, "1 and 2");

    state.SetComplexityN(state.range(0) * state.range(0));
}
BENCHMARK(VisualizeQuadCellData)->RangeMultiplier(2)->Range(2, 512)->Complexity();


//! benchmarks the VTU output
void VisualizeQuadWriteVtu(benchmark::State& state)
{
    using namespace NuTo::Visualize;
    TestCells s(state.range(0));

    Visualizer<VoronoiHandler> visu(s.cellInterfaces, VoronoiGeometryQuad(2));
    auto f = [&](const NuTo::CellData&, const NuTo::CellIpData& cipd) {
        auto coords = cipd.GlobalCoordinates();
        return Eigen::VectorXd::Constant(1, std::sin(10 * coords.x()) * std::cos(5 * coords.y()));
    };
    visu.DofValues(s.dof);
    visu.CellData(f, "fancy! - look at me!");

    for (auto _ : state)
        visu.WriteVtuFile("visuBenchmark" + std::to_string(state.range(0)) + ".vtu");

    state.SetComplexityN(state.range(0) * state.range(0));
}
BENCHMARK(VisualizeQuadWriteVtu)->RangeMultiplier(2)->Range(2, 512)->Complexity();


BENCHMARK_MAIN();
