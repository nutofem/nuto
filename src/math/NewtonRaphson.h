#pragma once

#include <Eigen/Sparse>
#include "base/Exception.h"
#include "math/LineSearch.h"

namespace NuTo
{
namespace NewtonRaphson
{

//! @brief custom exception for the newton algorithm
class NoConvergence : public NuTo::Exception
{
public:
    NoConvergence(const std::string& caller = "", const std::string& message = "")
        : Exception(caller, message)
    {
    }
};

//! @brief "Solver" for scalar values
struct DoubleSolver
{
    static double Solve(double dr, double r)
    {
        return r / dr;
    }
};

//! @brief Wrapper for Eigen::SparseMatrix
template <typename TSolver>
struct WrappedSolver
{
public:
    WrappedSolver(TSolver& solver)
        : mSolver(solver)
    {
    }

    Eigen::VectorXd Solve(const Eigen::SparseMatrix<double>& A, const Eigen::VectorXd& b)
    {
        mSolver.compute(A);
        return mSolver.solve(b);
    }

private:
    TSolver& mSolver;
};

template <typename TSolver>
WrappedSolver<TSolver> CreateWrappedSolver(TSolver& solver)
{
    return WrappedSolver<TSolver>(solver);
}


//! @brief problem definition
//! @tparam TR residual function
//! @tparam TDR derivative of the residual function
//! @tparam TNorm function that calculates the norm of the residual TTol = TNorm(TR)
//! @tparam TTol tolerance
//! @tparam TInfo Info function that takes (int, return type of TR, return type of TR)
template <typename TR, typename TDR, typename TNorm, typename TTol, typename TInfo>
struct Problem
{
    TR Residual;
    TDR Derivative;
    TNorm Norm;
    TTol mTolerance;
    TInfo Info;
};

//! @brief defines the problem, basically just to enable automatic template deduction. If you
//! create a Problem directly, you'll have to specify each template parameter. This methods avoids it.
template <typename TR, typename TDR, typename TNorm, typename TTol, typename TInfo = VoidInfo>
auto DefineProblem(TR residual, TDR derivative, TNorm norm, TTol tolerance, TInfo info = VoidInfo())
{
    return Problem<TR, TDR, TNorm, TTol, TInfo>({residual, derivative, norm, tolerance, info});
}

//! @brief solves the Problem using the newton raphson iteration with linesearch
//! @param problem type of the nonlinear problem
//! @param x0 of the initial value for the iteration
//! @param solver solver that provides a TX = solver.Solve(TNonlinearProblem::DR, TNonlinearProblem::R)
//! @param maxIterations default = 20
//! @param lineSearch line search algorithm, default = NoLineSearch, alternatively use NuTo::LineSearch()
//! @param numIterations optionally returns the number of iterations required
template <typename TNonlinearProblem, typename TX, typename TSolver, typename TLineSearchAlgorithm = NoLineSearch>
auto Solve(TNonlinearProblem&& problem, TX&& x0, TSolver&& solver, int maxIterations = 20,
           TLineSearchAlgorithm&& lineSearch = NoLineSearch(), int* numIterations = nullptr)
{
    auto x = x0;
    auto r = problem.Residual(x);

    int iteration = 0;
    problem.Info(iteration, x, r);

    if (problem.Norm(r) < problem.mTolerance)
        return x;

    while (iteration < maxIterations)
    {
        auto dr = problem.Derivative(x);
        auto dx = solver.Solve(dr, r);

        ++iteration;

        if (lineSearch(problem, &r, &x, dx))
        {
            if (numIterations)
                *numIterations = iteration;
            problem.Info(iteration, x, r);
            return x;
        }
        problem.Info(iteration, x, r);
    }
    if (numIterations)
        *numIterations = iteration;
    throw NoConvergence(__PRETTY_FUNCTION__, "No convergence after " + std::to_string(iteration) + " iterations.");
}
} /* NewtonRaphson */
} /* NuTo */
