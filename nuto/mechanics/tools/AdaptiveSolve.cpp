#include <iostream>
#include "nuto/base/Timer.h"
#include "nuto/base/Exception.h"
#include "nuto/math/NewtonRaphson.h"
#include "nuto/mechanics/tools/AdaptiveSolve.h"

#include <rang.hpp>

using namespace NuTo;

AdaptiveSolve::AdaptiveSolve(std::function<int(double)> doStepFunction, std::function<void(double)> postProcessFunction)
    : mDoStepFunction(doStepFunction)
    , mPostProcess(postProcessFunction)
{
}

void AdaptiveSolve::Solve(double tEnd, double tStart)
{
    Timer timer(__PRETTY_FUNCTION__, mPrintOutput);
    double t = tStart;
    if (mPrintOutput)
        std::cout << rang::fg::blue << "Starting adaptive solve from t = " << t << rang::style::reset << '\n';
    int iStep = 0;

    mPostProcess(t);

    while (true)
    {
        if (mPrintOutput)
            std::cout << "Step " << iStep << " at t = " << t << " with dt = " << dt << ".\n";

        try
        {
            int numIterations = mDoStepFunction(t + dt);

            if (mPrintOutput)
                std::cout << "Converence after " << numIterations << " iterations.\n";
            t += dt;
            iStep++;
            mPostProcess(t);

            if (t >= tEnd)
                break;

            if (numIterations < 3 && dt < dtMax)
            {
                dt *= increaseFactor;
                dt = std::min(dt, dtMax);
                if (mPrintOutput)
                    std::cout << rang::fg::green << rang::style::bold << "Increasing time step to " << dt
                              << rang::style::reset << '\n';
            }

            if (t + dt > tEnd)
                dt = tEnd - t;
        }
        catch (NewtonRaphson::NoConvergence& e)
        {
            dt *= decreaseFactor;
            if (mPrintOutput)
                std::cout << rang::fg::red << rang::style::bold << "Decreasing time step to " << dt
                          << rang::style::reset << '\n';

            if (dt < dtMin)
                throw Exception(__PRETTY_FUNCTION__,
                                "Time step smaller than prescribed minimal time step " + std::to_string(dtMin) + ".");
            continue; // without updating the global time
        }
    }
    if (mPrintOutput)
        std::cout << "Sucessfully reached time t = " << tEnd << "!\n";
}

void AdaptiveSolve::SetQuiet()
{
    mPrintOutput = false;
}
