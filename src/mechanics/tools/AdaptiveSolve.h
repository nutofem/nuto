#pragma once

#include <functional>

namespace NuTo
{
class AdaptiveSolve
{
public:
    //! ctor
    //! @param doStepFunction function that solves a problem for a given t and returns the number of iterations
    //! required. -1 indicates a failed solve
    //! @param postProcessFunction function that is called with the current time
    AdaptiveSolve(std::function<int(double)> doStepFunction,
                  std::function<void(double)> postProcessFunction = [](double) {});

    //! Perform an adaptive time integration until `tEnd` is reached
    //! @param tEnd end time
    void Solve(double tEnd);

    double dt = 0.1;
    double dtMax = 0.1;
    double dtMin = 1.e-6;

    double increaseFactor = 1.5;
    double decreaseFactor = 0.5;

private:
    std::function<int(double)> mDoStepFunction;
    std::function<void(double)> mPostProcess;
};
} /* NuTo */
