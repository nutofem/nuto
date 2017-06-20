#pragma once

#include "MathException.h"
#include "math/LineSearch.h"

namespace NuTo
{



//! @brief interface to solve a problem R(x) = 0 without the use of the derivative
//! @remark assumes that x and R are of the same type
template <typename TResidual, typename TNorm>
struct Residual
{
    typedef TNorm NormType;
    typedef TResidual ResidualType;

    virtual TNorm Norm(const TResidual&) = 0;
    virtual TResidual R(const TResidual&) = 0;
    virtual void Info(double iteration, const TResidual& x, const TResidual& r) const = 0;
};

//! @brief interface to solve a problem R(x) = 0 with the use of the derivative
template <typename TDerivative, typename TResidual, typename TNorm>
struct ResidualDerivative : Residual<TResidual, TNorm>
{
    typedef TDerivative DerivativeType;

    virtual TDerivative DR(const TResidual&) = 0;
};

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
    static double Solve (double dr, double r)
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

template <typename TNonlinearProblem, typename TX, typename TSolver>
auto Newton(TNonlinearProblem problem, TX x0, TSolver&& solver, int maxIterations = 20, int* numIterations = nullptr)
{
    auto x = x0;
    auto r = problem.ResidualFunction(x);

    int iteration = 0;
    while (iteration < maxIterations)
    {
        auto dr = problem.DerivativeFunction(x);
        x -= solver.Solve(dr, r);
        r = problem.ResidualFunction(x);
        ++iteration;

        if (problem.NormFunction(r) < problem.mTolerance)
        {
            if (numIterations)
                *numIterations = iteration;
            return x;
        }
    }
    throw NoConvergence(__PRETTY_FUNCTION__, "No convergence after " + std::to_string(iteration) + " iterations.");
}

//! @brief Newton-Raphson algrithm. https://en.wikipedia.org/wiki/Newton%27s_method
//! finds the root of TFunction::R(x) using the derivative TFunction::DR(x)
//! @tparam TFunction ... implementation of the NuTo::ResidualDerivative interface
//! @tparam TUseLinesearch ... option whether to use a line search optimization algorithm or not
template <typename TFunction, bool TUseLinesearch = false>
class NewtonRaphson
{
    using ResidualType = typename TFunction::ResidualType;
    using NormType = typename TFunction::NormType;
    using DerivativeType = typename TFunction::DerivativeType;

public:
    NewtonRaphson(NormType tolerance, int maxIterations = 20)
        : mLineSearch(tolerance)
        , mMaxIterations(maxIterations)
    {
    }

    //! @brief solves the nonlinear equation TFunction::R(x) = 0
    //! @param f ... defines function/derivative/norm
    //! @param x ... start value 
    //! @param solver ... solver that implements ResidualType Solve(DerivativeType, ResidualType)
    template <typename TSolver>
    ResidualType Solve(TFunction f, ResidualType x, TSolver&& solver, int* numIterations = nullptr) const
    {
        int iteration = 0;
        ResidualType r = f.R(x);
        f.Info(iteration, x, r);
        while (iteration < mMaxIterations)
        {
            ++iteration;
            DerivativeType dr = f.DR(x);
            ResidualType dx = solver.Solve(dr, r);

            bool acceptSolution = false;
            std::tie(acceptSolution, r, x) = mLineSearch.template Evaluate<TUseLinesearch>(x, dx, f);
            if (acceptSolution)
            {
                if (numIterations) 
                    *numIterations = iteration;
                return x;
            }

            f.Info(iteration, x, r);
        }
        throw NoConvergence(__PRETTY_FUNCTION__, "No convergence after " + std::to_string(iteration) + " iterations.");
    }

private:
    NuTo::LineSearch<TFunction> mLineSearch;
    int mMaxIterations;
};

} /* NuTo */
