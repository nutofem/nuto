#include <iomanip>
#include <iostream>
#include "optimize/Jacobi.h"
#include "base/Timer.h"

//! @brief ... Optimize routine - optimize displacement or error according to input
//! @brief ... Optimize residual if $|r|\inequal 0$
//! @brief ... Optimize displacements if |r|=0, forces have to be added (at place of residuals)
int NuTo::Jacobi::Optimize()
{
    std::vector<double>& v = GetParametersVec();
    std::vector<double>& f = mpCallbackHandlerGrid->GetRightHandSide();

    int returnValue = Optimize(v, f);

    return returnValue;
}

//! @brief ... Optimize routine - optimize solution vector (displacements or error)
//! @brief ... equation system: Kv=f or Ke=r
//! @brief ... $v \leftarrow v - \omega D^{-1} K v + \omega D^{-1} f
//! @param ... v - solution vector, f - right hand sight vector
//! @return ... optimization_return_attribute
int NuTo::Jacobi::Optimize(std::vector<double>& v, std::vector<double>& f)
{
    NuTo::Timer timer(__PRETTY_FUNCTION__, mShowTime);

    double rErrorNorm = 0.;
    int numGradientCalls(0), // number of gradient calls
            curIteration(0); // number of iterations

    eOptimizationReturnAttributes returnValue;


    std::vector<double> vNext(mNumParameters + 3, 0.);
    std::vector<double> rhs = f;

    int precision = 6;
    int width = 10;

    bool converged(false);

    int localMaxGradientCalls = (int)mNumParameters;
    //  int localMaxGradientCalls=2*mNumParameters;
    if (localMaxGradientCalls < mMaxGradientCalls)
        SetMaxGradientCalls(localMaxGradientCalls);

    mpCallbackHandlerGrid->Hessian(rhs);
    // weighting factor (normally at Hessian())
    for (size_t i = 0; i < mNumParameters; ++i)
        rhs[i] *= mOmega;

    while (!converged)
    {
        numGradientCalls++;
        if (numGradientCalls > mMaxGradientCalls)
        {
            converged = true;
            returnValue = eOptimizationReturnAttributes::MAXGRADIENTCALLS;
            break;
            // return MAXGRADIENTCALLS;
        }
        curIteration++;
        if (curIteration > mMaxIterations)
        {
            converged = true;
            returnValue = eOptimizationReturnAttributes::MAXITERATIONS;
            break;
            // return MAXGRADIENTCALLS;
        }
        // reset gNext to zero in gradient
        mpCallbackHandlerGrid->Gradient(v, vNext);
        // multiply with point diagonal preconditoner
        mpCallbackHandlerGrid->Hessian(vNext);
        rErrorNorm = 0.;
        for (size_t i = 0; i < mNumParameters; ++i)
        {
            // do not forget sign!!!
            vNext[i] *= -mOmega;
            vNext[i] += rhs[i];
            rErrorNorm += (vNext[i]) * (vNext[i]);
            v[i] += vNext[i];
        }

        if (rErrorNorm < mAccuracyGradient * mAccuracyGradient)
        {
            converged = true;
            returnValue = eOptimizationReturnAttributes::DELTAOBJECTIVEBETWEENCYCLES;
            break;
        }
    }

    mIsBuild = true;

    if (mVerboseLevel > 0)
    {
        std::cout << "[Jacobi::Optimize] Number of Iterations............. " << curIteration << std::endl;
        std::cout << "[Jacobi::Optimize] Active convergence criterion..... ";
        switch (returnValue)
        {
        case eOptimizationReturnAttributes::MAXGRADIENTCALLS:
            std::cout << "Maximum number of gradient calls reached." << std::endl;
            break;
        case eOptimizationReturnAttributes::MAXITERATIONS:
            std::cout << "Maximum number of iterations reached." << std::endl;
            break;
        case eOptimizationReturnAttributes::DELTAOBJECTIVEBETWEENCYCLES:
            std::cout << "Norm of error smaller than prescribed value." << std::endl;
            mObjective = sqrt(rErrorNorm);
            break;
        default:
            std::cout << "Unknown convergence criterion." << std::endl;
            break;
        }
        std::cout << std::endl;

        std::cout.precision(precision);
        std::cout << std::setw(width) << "[Jacobi::Optimize] parameters ";
        for (size_t count = 0; count < mNumParameters; count++)
        {
            std::cout << std::setw(width) << v[count] << "   ";
        }
        std::cout << std::endl;
    }
    return static_cast<int>(returnValue);
}

//! @brief ... Info routine that prints general information about the object (detail according to verbose level)
void NuTo::Jacobi::Info() const
{
    //    NuTo::Optimizer::InfoBase();
    std::cout << "AccuracyGradient" << mAccuracyGradient << std::endl;
    std::cout << "MinDeltaObjBetweenRestarts" << mMinDeltaObjBetweenRestarts << std::endl;
    std::cout << "MaxGradientCalls" << mMaxGradientCalls << std::endl;
    std::cout << "MaxHessianCalls" << mMaxHessianCalls << std::endl;
    std::cout << "MaxIterations" << mMaxIterations << std::endl;
    std::cout << "ShowSteps" << mShowSteps << std::endl;
}
