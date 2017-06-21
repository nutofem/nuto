#pragma once

namespace NuTo
{
//! @brief Performs the line search algorithm based on the results of a single newton iteration step
class LineSearch
{
public:
    //! @brief ctor
    //! @param maxNumLineSearchSteps ... maximal number of line search steps
    constexpr LineSearch(int maxNumLineSearchSteps = 6)
        : mMaxNumLineSearchStep(maxNumLineSearchSteps)
    {
    }

    //! @brief actual line search implementation
    //! @remark several return arguments are used to help inlining everything. A lot of time went into actual
    //! benchmarking. Dunno why, but the a previous version that avoids the return arguments and returns various values
    //! in a std::tuple was significantly slower (10%) in a benchmark for a scalar function
    //! @param problem ... class that implements NormFunction, ResidualFunction and mTolerance
    //! @param r ... residual, return argument, see remark
    //! @param x ... value of the argument x, return argument, see remark
    //! @param dx ... dx from the solver
    template <typename TProblem, typename TX>
    bool operator()(TProblem&& problem, TX* r, TX* x, TX dx) const
    {
        double alpha = 1.;
        int lineSearchStep = 0;
        const auto x0 = *x;
        const auto previousNorm = problem.NormFunction(problem.ResidualFunction(*x));
        while (lineSearchStep < mMaxNumLineSearchStep)
        {
            *x = x0 - alpha * dx;
            *r = problem.ResidualFunction(*x);
            const auto trialNorm = problem.NormFunction(*r);

            if (trialNorm < problem.mTolerance)
                return true;

            alpha *= 0.5;
            lineSearchStep++;

            if (trialNorm < (1. - alpha) * previousNorm)
                return false;
        }
        return false; // max steps reached
    }

private:
    int mMaxNumLineSearchStep;
};


//! @brief just a normal continuation of the newton scheme without using line search while keeping the interface of
//! NuTo::LineSearch
class NoLineSearch
{
public:
    template <typename TProblem, typename TX>
    bool operator()(TProblem&& problem, TX* r, TX* x, TX dx) const
    {
        *x -= dx;
        *r = problem.ResidualFunction(*x);
        return problem.NormFunction(*r) < problem.mTolerance;
    }
};
} /* NuTo */
