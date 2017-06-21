#include "Benchmark.h"
#include "math/NewtonRaphson.h"
#include <cmath>

constexpr double tolerance = 1.e-10;
constexpr double runtime = .1;

auto Problem()
{
    auto R = [](double x) { return x * x * x - x + 6; };
    auto DR = [](double x) { return 3. * x * x - 1; };
    auto Norm = [](double x) { return std::abs(x); };
    return NuTo::DefineNonlinearProblem(R, DR, Norm, tolerance); 
}

BENCHMARK(Newton, NuToFunction, runner)
{
    auto problem = Problem();
    //NuTo::NoLineSearch<decltype(problem), double> linesearch;
    //const auto linesearch = NuTo::NoLineSearch();
    while (runner.KeepRunningTime(runtime))
    {
        auto x = NuTo::Newton(problem, 0., NuTo::DoubleSolver(), 100);
        if (std::fabs(x + 2) > 1.e-10)
            throw;
    }
}

BENCHMARK(Newton, NuToFunctionLineSearch, runner)
{
    auto problem = Problem();
    while (runner.KeepRunningTime(runtime))
    {
        auto x = NuTo::Newton(problem, 0., NuTo::DoubleSolver(), 100, NuTo::LineSearch());
        if (std::fabs(x + 2) > 1.e-10)
            throw;
    }
}

BENCHMARK(Newton, hardcode, runner)
{
    while (runner.KeepRunningTime(runtime))
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
            if (std::fabs(r) < tolerance)
                break;
            if (i > 100)
                break;
        }
        if (std::fabs(x + 2) > 1.e-10)
            throw;
    }
}

