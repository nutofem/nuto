#pragma once

#include "MathException.h"

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
    virtual void Info(unsigned mIteration, const TResidual& x, const TResidual& r) const {};
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

template <typename TFunction>
class NewtonRaphson
{
    using ResidualType = typename TFunction::ResidualType;
    using NormType = typename TFunction::NormType;
    using DerivativeType = typename TFunction::DerivativeType;

public:
    NewtonRaphson(NormType tolerance, unsigned maxIterations = 20)
        : mTolerance(tolerance)
        , mMaxIterations(maxIterations)
    {
    }

    //! @brief solves the nonlinear equation R = 0
    template <typename TSolver>
    ResidualType Solve(TFunction f, ResidualType x, TSolver& solver) const
    {
        unsigned mIteration = 0;
        ResidualType r = f.R(x);
        f.Info(mIteration, x, r);
        while ((f.Norm(r) > mTolerance) & (mIteration < mMaxIterations))
        {
            DerivativeType dr = f.DR(x);
            x -= solver.Solve(dr, r);
            r = f.R(x);
            ++mIteration;
            f.Info(mIteration, x, r);
            if (f.Norm(r) < mTolerance)
                return x;
        }
        throw NoConvergence(__PRETTY_FUNCTION__, "No convergence after " + std::to_string(mIteration) + " iterations.");
    }

    // returns the number of iterations of the previous solve
    unsigned NumIterations() const
    {
        return mIteration;
    }

private:
    NormType mTolerance;
    unsigned mMaxIterations;
    unsigned mIteration = 0;
};

} /* NuTo */
