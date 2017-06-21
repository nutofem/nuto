#pragma once

#include "MathException.h"
#include "math/LineSearch.h"

namespace NuTo
{

class NoConvergence : public NuTo::MathException
{
public:
    NoConvergence(const std::string& caller, const std::string& message)
        : MathException(caller, message)
    {
    }
};

struct DoubleSolver
{
    static double Solve(double dr, double r)
    {
        return r / dr;
    }
};

template <typename TR, typename TDR, typename TNorm, typename TTol>
struct NonlinearProblem
{
    TR ResidualFunction;
    TDR DerivativeFunction;
    TNorm NormFunction;
    TTol mTolerance;
};

template <typename TR, typename TDR, typename TNorm, typename TTol>
auto DefineNonlinearProblem(TR residual, TDR derivative, TNorm norm, TTol tolerance)
{
    return NonlinearProblem<TR, TDR, TNorm, TTol>({residual, derivative, norm, tolerance});
}


template <typename TNonlinearProblem, typename TX, typename TSolver, typename TLineSearchAlgorithm = NoLineSearch>
auto Newton(TNonlinearProblem&& problem, TX&& x0, TSolver&& solver, int maxIterations = 20,
            TLineSearchAlgorithm&& lineSearch = NoLineSearch(), int* numIterations = nullptr)
{
    auto x = x0;
    auto r = problem.ResidualFunction(x);

    int iteration = 0;
    while (iteration < maxIterations)
    {
        const auto dr = problem.DerivativeFunction(x);
        const auto dx = solver.Solve(dr, r);

        ++iteration;
        if (lineSearch(problem, &r, &x, dx))
            return x;
    }
    throw NoConvergence(__PRETTY_FUNCTION__, "No convergence after " + std::to_string(iteration) + " iterations.");
}
} /* NuTo */
