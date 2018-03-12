#pragma once

#include "optimize/Optimizer.h"

#include <vector>

namespace NuTo
{

//! @author Joerg F. Unger, ISM
//! @date September 2009
//! @brief Conjugate gradient optimizer
class ConjugateGradientNonLinear : public Optimizer
{
public:
    ConjugateGradientNonLinear(unsigned int rNumParameters)
        : Optimizer(rNumParameters, (unsigned int)0, (unsigned int)0)
    {
        mAccuracyGradient = 1e-6;
        mMinDeltaObjBetweenRestarts = 1e-6;
        mMaxGradientCalls = INT_MAX, mMaxHessianCalls = INT_MAX, mMaxIterations = INT_MAX;
        mShowSteps = 1;
    }

    virtual ~ConjugateGradientNonLinear() = default;

    int Optimize() override;

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

    inline void SetShowSteps(int rShowSteps)
    {
        mShowSteps = rShowSteps;
    }

    //! @brief Info routine that prints general information about the object (detail according to verbose level)
    virtual void Info() const;

protected:
    void CalcScalingFactors(int& numHessianCalls, Eigen::MatrixXd& hessianOrig, Eigen::VectorXd& scaleFactorsInv);

    double mAccuracyGradient;
    double mMinDeltaObjBetweenRestarts;
    int mMaxGradientCalls;
    int mMaxHessianCalls;
    int mMaxIterations;
    int mShowSteps;
};
} // namespace NuTo
