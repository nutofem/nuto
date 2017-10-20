#pragma once

#include "optimize/Optimizer.h"

#include <vector>

namespace NuTo
{
//! @author Andrea Keszler, ISM
//! @date July 2010
//! @brief ... standard class for jacobi method
class Jacobi : public virtual Optimizer
{

public:
    Jacobi(unsigned int rNumParameters)
        : Optimizer(rNumParameters, (unsigned int)0, (unsigned int)0)
    {
        mAccuracyGradient = 1e-6;
        mMinDeltaObjBetweenRestarts = 1e-6;
        mOmega = 2. / 3.;
        mMaxGradientCalls = INT_MAX, mMaxHessianCalls = INT_MAX, mMaxIterations = INT_MAX;
        mShowSteps = 100;
        mNumParameters = rNumParameters;
    }

    virtual ~Jacobi() = default;


    int Optimize() override;

    int Optimize(std::vector<double>& v, std::vector<double>& f);

    inline void SetMaxGradientCalls(int rMaxGradientCalls)
    {
        mMaxGradientCalls = rMaxGradientCalls;
    }

    inline void SetMaxHessianCalls(int rMaxHessianCalls)
    {
        mMaxHessianCalls = rMaxHessianCalls;
    }

    inline void SetMaxIterations(int rMaxIterations)
    {
        mMaxIterations = rMaxIterations;
    }

    inline void SetAccuracyGradient(double rAccuracyGradient)
    {
        mAccuracyGradient = rAccuracyGradient;
    }

    inline void SetMinDeltaObjBetweenRestarts(double rMinDeltaObjBetweenRestarts)
    {
        mMinDeltaObjBetweenRestarts = rMinDeltaObjBetweenRestarts;
    }

    inline void SetScalingFactorOmega(double rOmega)
    {
        mOmega = rOmega;
    }

    inline void SetShowSteps(int rShowSteps)
    {
        mShowSteps = rShowSteps;
    }

    //! @brief ... Info routine that prints general information about the object (detail according to verbose level)
    virtual void Info() const;


protected:
    double mAccuracyGradient;
    double mMinDeltaObjBetweenRestarts;
    double mOmega; // scaling factor
    int mMaxGradientCalls;
    int mMaxHessianCalls;
    int mMaxIterations;
    int mShowSteps;
    size_t mNumParameters;
};
} // namespace NuTo
