#include <iostream>
#include "optimize/ConjugateGradientNonLinear.h"
#define machine_precision 1e-15
#define tol 1e-8 // sqrt machine_precision
#define goldenSect 0.38196601125
#define ZEPS 1e-10 // absolute tolerance
#define SHFT(a, b, c, d)                                                                                               \
    (a) = (b);                                                                                                         \
    (b) = (c);                                                                                                         \
    (c) = (d);
#define SIGN(a, b) ((b) >= 0.0 ? std::abs(a) : -std::abs(a))
#define golden_sec 0.38196601125

// @brief
// return reason
int NuTo::ConjugateGradientNonLinear::Optimize()
{
    double objectiveLastRestart, beta, initialAlpha(1e-2), v, w, fv, fw, e, d, u, prevAlpha, prevObjective,
            intermediateAlpha, intermediateObjective, alpha, normProjection, normPrevGradient, normGrad;

    int numFunctionCalls(0), // number of function calls
            numGradientCalls(0), // number of gradient calls
            numHessianCalls(0), // number of hessian calls
            curIteration(0), // number of iterations
            curCycle(0); // number of iterations without restart

    eOptimizationReturnAttributes returnValue;

    bool BrentsMethodConverged;

    Eigen::MatrixXd gradientOrig;
    Eigen::MatrixXd hessianOrig;
    Eigen::VectorXd prevParameters;
    Eigen::VectorXd gradientScaled;
    Eigen::VectorXd scaleFactorsInv(GetNumParameters());
    Eigen::VectorXd prevGradientScaled(GetNumParameters());
    Eigen::VectorXd deltaGradient(GetNumParameters());
    Eigen::VectorXd searchDirectionScaled(GetNumParameters());
    Eigen::VectorXd searchDirectionOrig(GetNumParameters());

    bool converged(false);
    double mAccuracyGradientScaled = mAccuracyGradient * sqrt(GetNumParameters());

    // check, if callback handler is set
    if (mpCallbackHandler == nullptr)
        throw Exception("[ConjugateGradientNonLinear::Optimize] Callback handler not set to determine "
                                "mObjective function and derivatives.");

    // calculate mObjective
    mObjective = mpCallbackHandler->Objective();
    objectiveLastRestart = mObjective;
    numFunctionCalls++;
    if (numFunctionCalls > mMaxFunctionCalls)
    {
        converged = true;
        returnValue = eOptimizationReturnAttributes::MAXHESSIANCALLS;
    }
    while (!converged)
    {
        if (mObjective < mMinObjective)
        {
            converged = true;
            returnValue = eOptimizationReturnAttributes::MINOBJECTIVE;
            break;
        }

        // calculate gradient
        mpCallbackHandler->Gradient(gradientOrig);
        numGradientCalls++;
        if (numGradientCalls > mMaxGradientCalls)
        {
            converged = true;
            returnValue = eOptimizationReturnAttributes::MAXGRADIENTCALLS;
            break;
            // return MAXGRADIENTCALLS;
        }

        // determine search direction
        if (curCycle % GetNumParameters() == 0)
        {
            // initialize search direction with steepest descent
            // calculate Hessian for scaling
            CalcScalingFactors(numHessianCalls, hessianOrig, scaleFactorsInv);
            if (numHessianCalls > mMaxHessianCalls)
            {
                converged = true;
                returnValue = eOptimizationReturnAttributes::MAXHESSIANCALLS;
                break;
            }
            gradientScaled = scaleFactorsInv.asDiagonal() * gradientOrig;

            normGrad = gradientScaled.norm();

            // printf("norm Grad %g\n",normGrad);
            if (normGrad < mAccuracyGradientScaled)
            {
                converged = true;
                returnValue = eOptimizationReturnAttributes::NORMGRADIENT;
                break;
            }

            searchDirectionScaled = -gradientScaled;
            searchDirectionOrig = scaleFactorsInv.asDiagonal() * searchDirectionScaled;
            if (mVerboseLevel > 5 && curCycle > 0)
                std::cout << "   Restart after " << curCycle << " cycles, DeltaObjective "
                          << objectiveLastRestart - mObjective << std::endl;
            objectiveLastRestart = mObjective;
            curCycle = 0;
        }
        else
        {
            // scale gradient
            gradientScaled = scaleFactorsInv.asDiagonal() * gradientOrig;

            normProjection = gradientScaled.dot(prevGradientScaled);
            normPrevGradient = prevGradientScaled.dot(prevGradientScaled);
            // Powell-Beale Restarts
            // printf("restart criterion %g>%g\n",std::abs(normProjection),0.2*normPrevGradient);
            if (std::abs(normProjection) > 0.2 * normPrevGradient)
            {
                // check for convergence
                if (std::abs(objectiveLastRestart - mObjective) / mObjective < mMinDeltaObjBetweenRestarts)
                {
                    converged = true;
                    returnValue = eOptimizationReturnAttributes::DELTAOBJECTIVEBETWEENCYCLES;
                    break;
                    // return DELTAOBJECTIVEBETWEENCYCLES;
                }
                else
                {
                    if (std::abs(objectiveLastRestart - mObjective) < mMinDeltaObjBetweenRestarts)
                    {
                        {
                            converged = true;
                            returnValue = eOptimizationReturnAttributes::DELTAOBJECTIVEBETWEENCYCLES;
                            break;
                        }
                    }
                }

                // restart with preconditioning
                // calculate hessian for preconditioning
                CalcScalingFactors(numHessianCalls, hessianOrig, scaleFactorsInv);
                if (numHessianCalls > mMaxHessianCalls)
                {
                    converged = true;
                    returnValue = eOptimizationReturnAttributes::MAXHESSIANCALLS;
                    break;
                }

                // scale gradient
                gradientScaled = scaleFactorsInv.asDiagonal() * gradientOrig;

                searchDirectionScaled = -gradientScaled;
                searchDirectionOrig = scaleFactorsInv.asDiagonal() * searchDirectionScaled;

                if (mVerboseLevel > 5)
                    std::cout << "   Restart after " << curCycle << " cycles, DeltaObjective "
                              << objectiveLastRestart - mObjective << std::endl;

                normGrad = gradientScaled.norm() * (double)GetNumParameters();

                objectiveLastRestart = mObjective;

                if (normGrad < mAccuracyGradientScaled)
                {
                    converged = true;
                    returnValue = eOptimizationReturnAttributes::NORMGRADIENT;
                    break;
                }
                curCycle = 0;
            }
            else
            {
                // update search direction according to Polak-Ribiere update formular
                deltaGradient = gradientScaled - prevGradientScaled;
                beta = deltaGradient.dot(gradientScaled) / normPrevGradient;

                if (beta < 0)
                {
                    // std::cout<< "Set beta ("<< beta <<") to zero " << std::endl;
                    beta = 0;
                }

                searchDirectionScaled *= beta;
                searchDirectionScaled -= gradientScaled;
                searchDirectionOrig = scaleFactorsInv.asDiagonal() * searchDirectionScaled;
                double normSearch = searchDirectionScaled.norm();

                if (normSearch < mAccuracyGradientScaled)
                {
                    // restart with preconditioning
                    // calculate hessian for preconditioning
                    CalcScalingFactors(numHessianCalls, hessianOrig, scaleFactorsInv);
                    if (numHessianCalls > mMaxHessianCalls)
                    {
                        converged = true;
                        returnValue = eOptimizationReturnAttributes::MAXHESSIANCALLS;
                        break;
                    }

                    // scale gradient
                    gradientScaled = scaleFactorsInv.asDiagonal() * gradientOrig;

                    searchDirectionScaled = -gradientScaled;
                    searchDirectionOrig = scaleFactorsInv.asDiagonal() * searchDirectionScaled;

                    if (mVerboseLevel > 3)
                        std::cout << "   Restart after " << curCycle << " cycles, DeltaObjective "
                                  << objectiveLastRestart - mObjective << std::endl;

                    objectiveLastRestart = mObjective;

                    normGrad = gradientScaled.norm();

                    if (normGrad < mAccuracyGradientScaled)
                    {
                        converged = true;
                        returnValue = eOptimizationReturnAttributes::NORMGRADIENT;
                        break;
                    }
                    curCycle = 0;
                }
            }
        }
        // store gradient and previous mObjective for next search
        prevGradientScaled = gradientScaled;

        /*
        int precision = 3;
        int width = 10;
        std::cout.precision(precision);
        std::cout << std::setw(width)<< "gradient scaled " ;
        for (int count=0; count<GetNumParameters(); count++)
        {
            std::cout << std::setw(width)<< gradientScaled(count) << "   " ;
        }
        std::cout << std::endl;
        std::cout << std::setw(width)<< "gradient orig " ;
        for (int count=0; count<GetNumParameters(); count++)
        {
            std::cout << std::setw(width)<< gradientOrig(count,0) << "   " ;
        }
        std::cout << std::endl;
        std::cout << std::setw(width)<< "searchDirectionScaled " ;
        for (int count=0; count<GetNumParameters(); count++)
        {
            std::cout << std::setw(width)<< searchDirectionScaled(count) << "   " ;
        }
        std::cout << std::endl;
        std::cout << std::setw(width)<< "searchDirectionOrig " ;
        for (int count=0; count<GetNumParameters(); count++)
        {
            std::cout << std::setw(width)<< searchDirectionOrig(count) << "   " ;
        }
        std::cout << std::endl;
        */

        if (mVerboseLevel > 1 && curIteration % mShowSteps == 0)
            std::cout << "Iteration " << curIteration << " with mObjective " << mObjective << " and norm grad scaled "
                      << gradientScaled.norm() / sqrt(GetNumParameters()) << std::endl;

        // increase iteration and curCycle
        curCycle++;
        curIteration++;

        if (curIteration > mMaxIterations)
        {
            converged = true;
            returnValue = eOptimizationReturnAttributes::MAXITERATIONS;
            break;
        }

        // store previous Parameters and mObjective for linesearch
        prevAlpha = 0;
        prevObjective = mObjective;
        intermediateAlpha = 0;
        intermediateObjective = mObjective;
        prevParameters = mvParameters;

        // perform line search
        alpha = initialAlpha / searchDirectionScaled.norm();

        mvParameters = prevParameters + alpha * searchDirectionOrig;

        // set new parameters in search direction
        mpCallbackHandler->SetParameters(mvParameters);

        // calculate mObjective
        mObjective = mpCallbackHandler->Objective();
        numFunctionCalls++;
        if (numFunctionCalls > mMaxFunctionCalls)
        {
            converged = true;
            returnValue = eOptimizationReturnAttributes::MAXFUNCTIONCALLS;
            break;
        }

        // find first interval with x1<x2<x3 and f(x1)>f(x2) und f(x2)<f(x3)
        int numFunctionCallsBefore(numFunctionCalls);
        if (mObjective < intermediateObjective)
        {
            if (mVerboseLevel > 5 && curIteration % mShowSteps == 0)
                std::cout << "   Increase search interval" << std::endl;
            while (mObjective < intermediateObjective)
            {
                // increase search interval until mObjective increases
                // printf("Increase interval with alphas %g %g %g and objectives %g %g
                // %g\n",prevAlpha,intermediateAlpha,alpha, prevObjective,intermediateObjective,mObjective);
                prevAlpha = intermediateAlpha;
                prevObjective = intermediateObjective;
                intermediateAlpha = alpha;
                intermediateObjective = mObjective;
                alpha += 1.5 * (intermediateAlpha - prevAlpha);
                mvParameters = prevParameters + alpha * searchDirectionOrig;
                mpCallbackHandler->SetParameters(mvParameters);
                mObjective = mpCallbackHandler->Objective();
                numFunctionCalls++;
                if (numFunctionCalls > mMaxFunctionCalls)
                {
                    converged = true;
                    returnValue = eOptimizationReturnAttributes::MAXFUNCTIONCALLS;
                    break;
                }
            }
        }
        else
        {
            if (mVerboseLevel > 5 && curIteration % mShowSteps == 0)
                std::cout << "   Decrease search interval" << std::endl;
            intermediateAlpha = alpha;
            intermediateObjective = mObjective;
            while (intermediateObjective >= prevObjective)
            {
                // decrease search interval until intermediate mObjective is smaller than previous mObjective
                // printf("Decrease interval with alphas %g %g %g and objectives %g %g
                // %g\n",prevAlpha,intermediateAlpha,alpha, prevObjective,intermediateObjective,mObjective);
                alpha = intermediateAlpha;
                mObjective = intermediateObjective;
                intermediateAlpha = 0.25 * alpha;
                mvParameters = prevParameters + intermediateAlpha * searchDirectionOrig;
                mpCallbackHandler->SetParameters(mvParameters);
                intermediateObjective = mpCallbackHandler->Objective();
                numFunctionCalls++;
                if (numFunctionCalls > mMaxFunctionCalls)
                {
                    converged = true;
                    returnValue = eOptimizationReturnAttributes::MAXFUNCTIONCALLS;
                    break;
                }
                if (std::abs(mObjective - intermediateObjective) < machine_precision)
                {
                    converged = true;
                    returnValue = eOptimizationReturnAttributes::REACHINGMACHINEPRECISION;
                    break;
                }
            }
        }
        if (mVerboseLevel > 3 && curIteration % mShowSteps == 0)
            std::cout << "   Number of iterations to determine interval " << numFunctionCalls - numFunctionCallsBefore
                      << std::endl;

        if (converged)
            break;
        // printf("Optimium is in between %g %g %g with alphas %g %g
        // %g\n",prevObjective,intermediateObjective,mObjective,prevAlpha,intermediateAlpha,alpha);

        // find optimum in interval given by prev_alpha intermediate_alpha alpha, with intermediateObjective smaller
        // than both interval limits
        v = intermediateAlpha;
        w = intermediateAlpha;
        fv = intermediateObjective;
        fw = intermediateObjective;
        d = goldenSect * (e = alpha - intermediateAlpha);
        u = intermediateAlpha + d;
        mvParameters = prevParameters + u * searchDirectionOrig;
        mpCallbackHandler->SetParameters(mvParameters);
        BrentsMethodConverged = false;
        numFunctionCallsBefore = numFunctionCalls;
        while (!BrentsMethodConverged)
        {
            mObjective = mpCallbackHandler->Objective();
            numFunctionCalls++;
            if (numFunctionCalls > mMaxFunctionCalls)
            {
                converged = true;
                returnValue = eOptimizationReturnAttributes::MAXFUNCTIONCALLS;
                break;
            }
            if (mObjective < intermediateObjective)
            {
                if (u >= intermediateAlpha)
                    prevAlpha = intermediateAlpha;
                else
                    alpha = intermediateAlpha;
                SHFT(v, w, intermediateAlpha, u);
                SHFT(fv, fw, intermediateObjective, mObjective);
            }
            else
            {
                if (u < intermediateAlpha)
                    prevAlpha = u;
                else
                    alpha = u;
                if (mObjective <= fw || w == intermediateAlpha)
                {
                    v = w;
                    w = u;
                    fv = fw;
                    fw = mObjective;
                }
                else
                {
                    if (mObjective <= fv || v == intermediateAlpha || v == w)
                    {
                        v = u;
                        fv = mObjective;
                    }
                }
            }

            // Done with housekeeping. Back for another iteration
            double tol1, tol2;
            double xm = 0.5 * (alpha + prevAlpha);
            tol2 = 2.0 * (tol1 = tol * std::abs(intermediateAlpha + ZEPS));
            if (std::abs(intermediateAlpha - xm) <= tol2 - 0.5 * (alpha - prevAlpha))
            {
                // line search converged
                mObjective = intermediateObjective;
                mvParameters = prevParameters + intermediateAlpha * searchDirectionOrig;
                initialAlpha = searchDirectionScaled.norm() * intermediateAlpha;
                BrentsMethodConverged = true;
                break;
            }
            else
            {
                if (std::abs(e) > tol1)
                {
                    // Construct a trial parabolic t.
                    double r = (intermediateAlpha - w) * (intermediateObjective - fv);
                    double q = (intermediateAlpha - v) * (intermediateObjective - fw);
                    double p = (intermediateAlpha - v) * q - (intermediateAlpha - w) * r;
                    q = 2.0 * (q - r);
                    if (q > 0.0)
                    {
                        p = -p;
                    }
                    q = std::abs(q);
                    double etemp = e;
                    e = d;
                    if (std::abs(p) >= std::abs(0.5 * q * etemp) || p <= q * (prevAlpha - intermediateAlpha) ||
                        p >= q * (alpha - intermediateAlpha))
                    {
                        d = golden_sec *
                            (e = (intermediateAlpha >= xm ? prevAlpha - intermediateAlpha : alpha - intermediateAlpha));
                    }
                    // The above conditions determine the acceptability of the parabolic t. Here we
                    // take the golden section step into the larger of the two segments.
                    else
                    {
                        d = p / q; // Take the parabolic step.
                        u = intermediateAlpha + d;
                        if (u - prevAlpha < tol2 || alpha - u < tol2)
                            d = SIGN(tol1, xm - intermediateAlpha);
                    }
                }
                else
                {
                    d = golden_sec *
                        (e = (intermediateAlpha >= xm ? prevAlpha - intermediateAlpha : alpha - intermediateAlpha));
                }
                u = (std::abs(d) >= tol1 ? intermediateAlpha + d : intermediateAlpha + SIGN(tol1, d));
                mvParameters = prevParameters + u * searchDirectionOrig;
                mpCallbackHandler->SetParameters(mvParameters);
            }
        }
        if (mVerboseLevel > 3 && curIteration % mShowSteps == 0)
            std::cout << "   Number of function calls for linesearch " << numFunctionCalls - numFunctionCallsBefore
                      << std::endl;
        if (converged)
            break;
    }
    mIsBuild = true;
    if (mVerboseLevel > 0)
    {
        std::cout << "Number of Function Calls......... " << numFunctionCalls << std::endl;
        std::cout << "Number of Gradient Calls......... " << numGradientCalls << std::endl;
        std::cout << "Number of Hessian Calls.......... " << numHessianCalls << std::endl;
        std::cout << "Number of Iterations............. " << curIteration << std::endl;
        std::cout << "Final mObjective function......... " << mObjective << std::endl;
        std::cout << "Norm of preconditioned gradient.. " << gradientScaled.norm() / sqrt(GetNumParameters())
                  << std::endl;
        std::cout << "Active convergence criterion..... ";
        switch (returnValue)
        {
        case eOptimizationReturnAttributes::MAXFUNCTIONCALLS:
            std::cout << "Maximum number of function calls reached." << std::endl;
            break;
        case eOptimizationReturnAttributes::MAXGRADIENTCALLS:
            std::cout << "Maximum number of gradient calls reached." << std::endl;
            break;
        case eOptimizationReturnAttributes::MAXHESSIANCALLS:
            std::cout << "Maximum number of hessian calls reached." << std::endl;
            break;
        case eOptimizationReturnAttributes::MAXITERATIONS:
            std::cout << "Maximum number of iterations reached." << std::endl;
            break;
        case eOptimizationReturnAttributes::NORMGRADIENT:
            std::cout << "Norm of preconditioned gradient smaller than prescribed value." << std::endl;
            break;
        case eOptimizationReturnAttributes::MINOBJECTIVE:
            std::cout << "Objective smaller than prescribed value." << std::endl;
            break;
        case eOptimizationReturnAttributes::DELTAOBJECTIVEBETWEENCYCLES:
            std::cout << "Delta mObjective between two restarts of the conjugate gradient method is smaller than "
                         "prescribed minimum."
                      << std::endl;
            break;
        case eOptimizationReturnAttributes::REACHINGMACHINEPRECISION:
            std::cout << "The machine precision is reached within the initial phase of the linesearch." << std::endl;
            break;
        default:
            std::cout << "Unknown convergence criterion." << std::endl;
        }
        std::cout << std::endl;
    }
    return static_cast<int>(returnValue);
}

void NuTo::ConjugateGradientNonLinear::CalcScalingFactors(int& numHessianCalls, Eigen::MatrixXd& hessianOrig,
                                                          Eigen::VectorXd& scaleFactorsInv)
{
    // calculate hessian for preconditioning
    mpCallbackHandler->Hessian(hessianOrig);
    // hessianOrig.Info();
    numHessianCalls++;

    // determine scale factors from the diagonal entries of the hessian
    for (int count = 0; count < GetNumParameters(); count++)
    {
        if (hessianOrig(count, count) > 1)
        {
            scaleFactorsInv(count) = 1. / hessianOrig(count, count);
        }
        else
        {
            scaleFactorsInv(count) = 1.;
        }
    }
}


//! @brief ... Info routine that prints general information about the object (detail according to verbose level)
void NuTo::ConjugateGradientNonLinear::Info() const
{
    NuTo::Optimizer::InfoBase();
    std::cout << "AccuracyGradient" << mAccuracyGradient << std::endl;
    std::cout << "MinDeltaObjBetweenRestarts" << mMinDeltaObjBetweenRestarts << std::endl;
    std::cout << "MaxGradientCalls" << mMaxGradientCalls << std::endl;
    std::cout << "MaxHessianCalls" << mMaxHessianCalls << std::endl;
    std::cout << "MaxIterations" << mMaxIterations << std::endl;
    std::cout << "ShowSteps" << mShowSteps << std::endl;
    std::cout << "AccuracyGradient" << mAccuracyGradient << std::endl;
}
