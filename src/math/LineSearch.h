#pragma once

namespace NuTo
{
namespace NewtonRaphson
{

//! @brief takes any argument and does nothing...
struct VoidInfo
{
    template <typename... TAny>
    void operator()(TAny&&...) const
    {
    }
};

//! @brief Performs the line search algorithm based on the results of a single newton iteration step
template <typename TInfo>
class LineSearchImplementation
{
public:
    //! @brief ctor
    //! @param maxNumLineSearchSteps maximal number of line search steps
    constexpr LineSearchImplementation(TInfo info, int maxNumLineSearchSteps)
        : mInfo(info)
        , mMaxNumLineSearchStep(maxNumLineSearchSteps)
    {
    }

    //! @brief actual line search implementation
    //! @remark several return arguments are used to help inlining everything. A lot of time went into actual
    //! benchmarking. Dunno why, but the a previous version that avoids the return arguments and returns various values
    //! in a std::tuple was significantly slower (10%) in a benchmark for a scalar function
    //! @param problem class that implements NormFunction, ResidualFunction and mTolerance
    //! @param r residual, return argument, see remark
    //! @param x value of the argument x, return argument, see remark
    //! @param dx dx from the solver
    template <typename TProblem, typename TX>
    bool operator()(TProblem&& problem, TX* r, TX* x, TX dx) const
    {
        double alpha = 1.;
        int lineSearchStep = 0;
        const auto x0 = *x;
        const auto previousNorm = problem.Norm(problem.Residual(*x));
        while (lineSearchStep < mMaxNumLineSearchStep)
        {
            *x = x0 - alpha * dx;
            *r = problem.Residual(*x);
            const auto trialNorm = problem.Norm(*r);

            mInfo(lineSearchStep, alpha, trialNorm);

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
    TInfo mInfo;
    int mMaxNumLineSearchStep;
};

//! @brief convienient instantiation of the NuTo::LineSearchImplementation with template deduction
template <typename TInfo = VoidInfo>
LineSearchImplementation<TInfo> LineSearch(TInfo info = VoidInfo(), int mMaxNumLineSearchStep = 6)
{
    return LineSearchImplementation<TInfo>(info, mMaxNumLineSearchStep);
}

//! @brief just a normal continuation of the newton scheme without using line search while keeping the interface of
//! NuTo::LineSearch
class NoLineSearch
{
public:
    template <typename TProblem, typename TX>
    bool operator()(TProblem&& problem, TX* r, TX* x, TX dx) const
    {
        *x -= dx;
        *r = problem.Residual(*x);
        return problem.Norm(*r) < problem.mTolerance;
    }
};

} /* NewtonRaphson */
} /* NuTo */
