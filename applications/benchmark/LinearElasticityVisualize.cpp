#include "Benchmark.h"
#include "LinearElasticBenchmarkStructure.h"

void ApplyRandomDisplacements(NuTo::Structure& s)
{
    auto allNodeIds = s.GroupGetMemberIds(s.GroupGetNodesTotal());
    for (auto id : allNodeIds)
        s.NodeSetDisplacements(id, Eigen::Vector3d::Random());
}

BENCHMARK(LinearElasticity, VisualizeBinary, runner)
{
    std::vector<int> numElements{10, 10, 100}; // i.e. 10k Elements
    NuTo::Benchmark::LinearElasticBenchmarkStructure s(numElements);
    ApplyRandomDisplacements(s.GetStructure());
    s.SetupVisualization();
    while(runner.KeepRunningIterations(10))
    {
        s.GetStructure().ExportVtkDataFileElements("LinearElasticVisualize.vtu");
    }
}

BENCHMARK(LinearElasticity, VisualizeAscii, runner)
{
    std::vector<int> numElements{10, 10, 100}; // i.e. 10k Elements
    NuTo::Benchmark::LinearElasticBenchmarkStructure s(numElements);
    ApplyRandomDisplacements(s.GetStructure());
    s.SetupVisualization();
    while(runner.KeepRunningIterations(10))
    {
        s.GetStructure().ExportVtkDataFileElements("LinearElasticVisualizeAscii.vtu", false);
    }
}
