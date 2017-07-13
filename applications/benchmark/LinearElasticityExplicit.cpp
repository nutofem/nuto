#include "Benchmark.h"
#include "LinearElasticBenchmarkStructure.h"
#include "mechanics/timeIntegration/RungeKutta4.h"

BENCHMARK(LinearElasticity, CompleteRun, runner)
{
    std::vector<int> numElements{10, 10, 100}; // i.e. 10k Elements
    NuTo::Benchmark::LinearElasticBenchmarkStructure s(numElements);
    NuTo::RungeKutta4 rk4(&s.GetStructure());
    rk4.SetTimeStep(1);
    rk4.SetResultDirectory("LinearElasticityExplicitResults", true);

    while (runner.KeepRunningIterations(1))
    {
        rk4.Solve(10);
    }
}
