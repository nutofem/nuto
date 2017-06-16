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

struct DoubleSolver
{
    double Solve(double dr, double r) const
    {
        return r / dr;
    }
};

constexpr double tolerance = 1.e-10;
constexpr double runtime = 0.1;

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
    DoubleSolver solver;
    while (runner.KeepRunningTime(runtime))
    {
        newton.Solve(f, 0, solver);
    }
}

BENCHMARK(Newton, NuToLineSearch, runner)
{
    F f;
    NuTo::NewtonRaphson<F, true> newton(tolerance, 100);
    DoubleSolver solver;
    while (runner.KeepRunningTime(runtime))
    {
        newton.Solve(f, 0, solver);
    }
}
