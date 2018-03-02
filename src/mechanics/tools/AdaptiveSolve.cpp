#include <iostream>
#include "base/Timer.h"
#include "base/Exception.h"
#include "math/NewtonRaphson.h"
#include "mechanics/tools/AdaptiveSolve.h"

using namespace NuTo;

AdaptiveSolve::AdaptiveSolve(std::function<int(double)> doStepFunction, std::function<void(double)> postProcessFunction)
    : mDoStepFunction(doStepFunction)
    , mPostProcess(postProcessFunction)
{
}

void AdaptiveSolve::Solve(double tEnd)
{
    Timer timer(__PRETTY_FUNCTION__, true);
    double t = 0;
    int iStep = 0;

    mPostProcess(t);

    while (true)
    {
        std::cout << "Step " << iStep << " at t = " << t << " with dt = " << dt << ".\n";

        try
        {
            int numIterations = mDoStepFunction(t + dt);

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
                std::cout << "--> Increasing time step to " << dt << ".\n";
            }

            if (t + dt > tEnd)
                dt = tEnd - t;
        }
        catch (...)
        {
            // std::cout << e.what() << '\n';
            dt *= decreaseFactor;
            std::cout << "--> Decreasing time step to " << dt << ".\n";

            if (dt < dtMin)
                throw Exception(__PRETTY_FUNCTION__,
                                "Time step smaller than prescribed minimal time step " + std::to_string(dtMin) + ".");
            continue; // without updating the global time
        }
    }
    std::cout << "Sucessfully reached time t = " << tEnd << "!\n";
}
