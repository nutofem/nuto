#include "Benchmark.h"
#include "math/NewtonRaphson.h"
#include <cmath>


struct F : NuTo::ResidualDerivative<double, double, double>
{
    double Norm(const double& x) override
    {
        return std::fabs(x);
    }
    double R(const double& x) override
    {
        return x * x * x - x + 6;
    }
    double DR(const double& x) override
    {
        return 3. * x * x - 1;
    }
    void Info(double, const double&, const double&) const override
    {
    }
};


constexpr double tolerance = 1.e-10;
constexpr double runtime = 0.1;

BENCHMARK(Newton, NuToFunction, runner)
{
    auto R = [](double x) { return x * x * x - x + 6; };
    auto DR = [](double x) { return 3. * x * x - 1; };
    auto Norm = [](double x) { return std::abs(x); };

    auto problem = NuTo::DefineNonlinearProblem(R, DR, Norm, tolerance);
    while (runner.KeepRunningTime(runtime))
    {
        auto x = NuTo::Newton(problem, 0., NuTo::DoubleSolver(), 100);
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


BENCHMARK(Newton, NuTo, runner)
{
    F f;
    NuTo::NewtonRaphson<F> newton(tolerance, 100);
    while (runner.KeepRunningTime(runtime))
    {
        newton.Solve(f, 0, NuTo::DoubleSolver());
    }
}

BENCHMARK(Newton, NuToLineSearch, runner)
{
    F f;
    NuTo::NewtonRaphson<F, true> newton(tolerance, 100);
    while (runner.KeepRunningTime(runtime))
    {
        newton.Solve(f, 0, NuTo::DoubleSolver());
    }
}
