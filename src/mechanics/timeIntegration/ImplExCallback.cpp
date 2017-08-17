//
// Created by Thomas Titscher on 1/17/17.
//

#include <cmath>
#include "mechanics/timeIntegration/ImplExCallback.h"

bool NuTo::ImplExCallback::AcceptSolution(double rMaxError, double rThreshold)
{
    bool extrapolationIsOK = rMaxError < rThreshold;
    if (not extrapolationIsOK && not mForceAcceptOfNextSolution)
    {
        mForceAcceptOfNextSolution = true;
        return false;
    }
    // accept solution, either the extrapolation is OK or this solution is forced to be accepted (somehow continue the
    // time stepping)
    mForceAcceptOfNextSolution = false;
    return true;
}

double NuTo::ImplExCallback::GetNewTimeStep(double rMaxError, double rThreshold, double mOldTimeStep)
{
    double newTimeStep = (rThreshold / rMaxError) * mOldTimeStep;
    return newTimeStep;
}

double NuTo::ImplExCallback::GetError(double, double k_n, double k_nn) const
{
    return std::abs(k_n - k_nn);
}
