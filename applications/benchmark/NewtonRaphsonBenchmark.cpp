#include "Benchmark.h"
#include "math/NewtonRaphson.h"
#include <cmath>

constexpr double tolerance = 1.e-10;
constexpr double runtime = .1;

auto ValidProblem()
{
    auto R = [](double x) { return x * x * x - x + 6; };
    auto DR = [](double x) { return 3. * x * x - 1; };
    auto Norm = [](double x) { return std::abs(x); };
    return NuTo::NewtonRaphson::DefineProblem(R, DR, Norm, tolerance);
}

auto InfoProblem()
{
    auto R = [](double x) { return x * x * x - x + 6; };
    auto DR = [](double x) { return 3. * x * x - 1; };
    auto Norm = [](double x) { return std::abs(x); };
    auto Info = [](int i, double x, double r) { std::cout << "NewtonStep: " << i << '\t' << x << '\t' << r << '\n'; };
    return NuTo::NewtonRaphson::DefineProblem(R, DR, Norm, tolerance, Info);
}

BENCHMARK(Newton, NuToFunction, runner)
{
    auto problem = ValidProblem();
    while (runner.KeepRunningTime(runtime))
    {
        auto x = NuTo::NewtonRaphson::Solve(problem, 0., NuTo::NewtonRaphson::DoubleSolver(), 100);
        if (std::fabs(x + 2) > 1.e-10)
            throw;
    }
}

BENCHMARK(Newton, NuToFunctionWithInfoToCout, runner)
{
    auto problem = InfoProblem();
    while (runner.KeepRunningTime(runtime))
    {
        auto x = NuTo::NewtonRaphson::Solve(problem, 0., NuTo::NewtonRaphson::DoubleSolver(), 100);
        if (std::fabs(x + 2) > 1.e-10)
            throw;
    }
}


BENCHMARK(Newton, NuToFunctionLineSearch, runner)
{
    auto problem = ValidProblem();
    while (runner.KeepRunningTime(runtime))
    {
        auto x = NuTo::NewtonRaphson::Solve(problem, 0., NuTo::NewtonRaphson::DoubleSolver(), 100,
                                            NuTo::NewtonRaphson::LineSearch());
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
