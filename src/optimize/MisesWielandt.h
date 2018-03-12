#pragma once

#include "optimize/Optimizer.h"
#include "base/Exception.h"

namespace NuTo
{
//! @author Andrea Keszler, ISM
//! @date July 2010
//! @brief standard class for vonMises-Wielandt method
class MisesWielandt : public virtual Optimizer
{

public:
    MisesWielandt(unsigned int rNumParameters)
        : Optimizer(rNumParameters, (unsigned int)0, (unsigned int)0)
    {
        mAccuracyGradient = 1e-6;
        mMinDeltaObjBetweenRestarts = 1e-6;
        mMaxGradientCalls = INT_MAX, mMaxHessianCalls = INT_MAX, mMaxIterations = INT_MAX;
        mShowSteps = 100;
        mUseDiagHessian = true;
        mNumParameters = rNumParameters;
        mObjectiveType = CONDITION_NUMBER_OF_MATRIX;
    }

    virtual ~MisesWielandt() = default;

    enum eObjectiveType
    {
        MAX_EIGENVALUE_OF_PRECOND_MATRIX,
        MAX_EIGENVALUE_OF_MATRIX,
        SPECTRAL_RADIUS_OF_MATRIX,
        SPECTRAL_RADIUS_OF_PRECOND_MATRIX,
        CONDITION_NUMBER_OF_MATRIX,
        CONDITION_NUMBER_OF_PRECOND_MATRIX,
    };

    int Optimize() override;


    void SetObjectiveType(std::string rObjectiveType)
    {
        std::string upperCaseObjectiveType;
        std::transform(rObjectiveType.begin(), rObjectiveType.end(), std::back_inserter(upperCaseObjectiveType),
                       (int (*)(int))toupper);

        if (upperCaseObjectiveType == "MAX_EIGENVALUE_OF_PRECOND_MATRIX")
            mObjectiveType = MAX_EIGENVALUE_OF_PRECOND_MATRIX;
        else if (upperCaseObjectiveType == "MAX_EIGENVALUE_OF_MATRIX")
            mObjectiveType = MAX_EIGENVALUE_OF_MATRIX;
        else if (upperCaseObjectiveType == "SPECTRAL_RADIUS_OF_MATRIX")
            mObjectiveType = SPECTRAL_RADIUS_OF_MATRIX;
        else if (upperCaseObjectiveType == "SPECTRAL_RADIUS_OF_PRECOND_MATRIX")
            mObjectiveType = SPECTRAL_RADIUS_OF_PRECOND_MATRIX;
        else if (upperCaseObjectiveType == "CONDITION_NUMBER_OF_MATRIX")
            mObjectiveType = CONDITION_NUMBER_OF_MATRIX;
        else if (upperCaseObjectiveType == "CONDITION_NUMBER_OF_PRECOND_MATRIX")
            mObjectiveType = CONDITION_NUMBER_OF_PRECOND_MATRIX;
        else
        {
            throw Exception("[NuTo::MisesWielandt::SetObjectiveType] ObjectiveType " + upperCaseObjectiveType +
                            " does not exist.");
        }
    }

    inline void SetObjectiveType(eObjectiveType rObjectiveType)
    {
        mObjectiveType = rObjectiveType;
    }
    inline eObjectiveType GetObjectiveType()
    {
        return mObjectiveType;
    }

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
    double mAccuracyGradient;
    double mMinDeltaObjBetweenRestarts;
    int mMaxGradientCalls;
    int mMaxHessianCalls;
    int mMaxIterations;
    int mShowSteps;
    bool mUseDiagHessian;
    eObjectiveType mObjectiveType;
    size_t mNumParameters;
};
} // namespace NuTo
