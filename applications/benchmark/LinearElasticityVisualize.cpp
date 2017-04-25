#include "Benchmark.h"
#include "LinearElasticBenchmarkStructure.h"

BENCHMARK(LinearElasticity, Visualize, runner)
{
    std::vector<int> numElements{10, 10, 100}; // i.e. 10k Elements
    NuTo::Benchmark::LinearElasticBenchmarkStructure s(numElements);
    s.SetupVisualization();
    while(runner.KeepRunningIterations(10))
    {
        s.GetStructure().ExportVtkDataFileElements("LinearElasticVisualize.vtu"); 
    }
}
