#include "Benchmark.h"
#include "LinearElasticBenchmarkStructure.h"

BENCHMARK(LinearElasticity, CompleteRun, runner)
{
    std::vector<int> numElements{10, 10, 100}; // i.e. 100k Elements
    while(runner.KeepRunningIterations(1))
    {
        NuTo::Benchmark::LinearElasticBenchmarkStructure s(numElements);
        s.SetupVisualization();
        s.SolveWithNewmark(10, "LinearElasticityResults");
    }
}
