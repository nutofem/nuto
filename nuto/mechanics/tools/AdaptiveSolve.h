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

    //! Perform an adaptive time integration from `t = tStart` to `t = tEnd`.
    //! One time step starts at time t and tries to apply the increment dt. If the step finishes without an exception,
    //! it is finalized by updating the time to t+dt calling the postprocess function. Only NoConvergence exception will
    //! be caught. They cause a decrease in the time step.
    //! @param tEnd end time
    //! @param tStart start time
    //! @remark tStart comes second (may be unintuitive) to default it to zero.
    void Solve(double tEnd, double tStart = 0);

    void SetQuiet();

    double dt = 0.1;
    double dtMax = 0.1;
    double dtMin = 1.e-6;

    double increaseFactor = 1.5;
    double decreaseFactor = 0.5;

private:
    std::function<int(double)> mDoStepFunction;
    std::function<void(double)> mPostProcess;
    bool mPrintOutput = true;
};
} /* NuTo */
