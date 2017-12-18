#include <cmath>
#include <benchmark/benchmark.h>

#include "math/NewtonRaphson.h"

/*
 * Shows that the performance of the NuTo::NewtonRaphson algorithm is very close to the performace of a handwritten
 * hardcode algorithm for a cubic scalar equation.
 */

constexpr double tolerance = 1.e-10;

auto ValidProblem()
{
    auto R = [](double x) { return x * x * x - x + 6; };
    auto DR = [](double x) { return 3. * x * x - 1; };
    auto Norm = [](double x) { return std::abs(x); };
    return NuTo::NewtonRaphson::DefineProblem(R, DR, Norm, tolerance);
}

void Check(double x)
{
    if (std::fabs(x + 2) > 1.e-10)
    {
        throw NuTo::NewtonRaphson::NoConvergence();
    }
}

static void NuToAlgorithm(benchmark::State& state)
{
    auto problem = ValidProblem();
    for (auto _ : state)
    {
        auto x = NuTo::NewtonRaphson::Solve(problem, 0., NuTo::NewtonRaphson::DoubleSolver(), 100);
        Check(x);
    }
}
BENCHMARK(NuToAlgorithm);

static void NuToWithLineSearch(benchmark::State& state)
{
    auto problem = ValidProblem();
    for (auto _ : state)
    {
        auto x = NuTo::NewtonRaphson::Solve(problem, 0., NuTo::NewtonRaphson::DoubleSolver(), 100,
                                            NuTo::NewtonRaphson::LineSearch());
        Check(x);
    }
}
BENCHMARK(NuToWithLineSearch);

static void Hardcode(benchmark::State& state)
{
    for (auto _ : state)
    {
        double x = 0;
        double r = x * x * x - x + 6;
        int i = 0;
        while (true)
        {
            double dr = 3. * x * x - 1;
            x -= r / dr;
            r = x * x * x - x + 6;
            ++i;
            if (std::fabs(r) < tolerance || i > 100)
            {
                break;
            }
        }
        Check(x);
    }
}
BENCHMARK(Hardcode);

BENCHMARK_MAIN();
