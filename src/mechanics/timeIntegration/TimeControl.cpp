#include "TimeControl.h"

#include "base/Exception.h"

void NuTo::TimeControl::Proceed()
{
    if (mTimeStep <= 0.0)
        throw Exception(__PRETTY_FUNCTION__, "Current timestep is 0 or negative!");
    mPreviousTime = mCurrentTime;
    mCurrentTime += mTimeStep;
}

void NuTo::TimeControl::AdjustTimestep(int iterations, int maxIterations, bool converged)
{


    mTimeStep = mTimeStepFunction(*this, iterations, maxIterations, converged);

    if (mTimeStep > mMaxTimeStep)
        mTimeStep = mMaxTimeStep;

    if (not converged && mPreviousTime != mCurrentTime)
        throw Exception(__PRETTY_FUNCTION__, "No convergence with the current maximum number of "
                                             "iterations, either use automatic time stepping, "
                                             "reduce the time step or the minimal line search cut "
                                             "back factor. In case you provided a custom timestepping function "
                                             "and intend to use some kind of automatic timestepping, call "
                                             "the RestorePreviosTime() function of the time control before reducing "
                                             "the timestep.");
}


void NuTo::TimeControl::SetTimeStep(double timeStep)
{
    if (timeStep <= 0)
        throw Exception(__PRETTY_FUNCTION__, "Timestep must be larger than 0!");

    mTimeStep = timeStep;
}

void NuTo::TimeControl::SetTimeFinal(double timeFinal)
{
    if (timeFinal <= mCurrentTime)
    {
        throw Exception(__PRETTY_FUNCTION__, "Final time must be larger than current time!");
    }
    mTimeFinal = timeFinal;
}


void NuTo::TimeControl::SetTimeStepFunction(std::function<double(TimeControl&, int, int, bool)> TimeStepFunction)
{
    mTimeStepFunction = TimeStepFunction;
}

void NuTo::TimeControl::UseDefaultAutomaticTimestepping()
{

    SetTimeStepFunction(DefaultAutomaticTimestepFunction);
}

void NuTo::TimeControl::UseEquidistantTimestepping()
{
    SetTimeStepFunction([](TimeControl& rTimeControl, int iterations, int maxIterations, bool converged) -> double {
        return rTimeControl.GetTimeStep();
    });
}

void NuTo::TimeControl::SetMaxTimeStep(double rMaxTimeStep)
{
    if (rMaxTimeStep <= 0.0)
        throw Exception(__PRETTY_FUNCTION__, "Maximal timestep must be a positive number!");
    if (rMaxTimeStep < mMinTimeStep)
        throw Exception(__PRETTY_FUNCTION__, "Maximal timestep must be bigger than minimal Timestep!");
    mMaxTimeStep = rMaxTimeStep;
}

void NuTo::TimeControl::SetMinTimeStep(double rMinTimeStep)
{
    if (rMinTimeStep > mMaxTimeStep)
        throw Exception(__PRETTY_FUNCTION__, "Minimal timestep must be smaller than maximal Timestep!");
    mMinTimeStep = rMinTimeStep;
}

double NuTo::TimeControl::DefaultAutomaticTimestepFunction(TimeControl& rTimeControl, int iterations, int maxIterations,
                                                           bool converged)
{
    if (converged)
    {
        if (iterations < 0.25 * maxIterations)
            return rTimeControl.GetTimeStep() * 1.5;
        else
            return rTimeControl.GetTimeStep();
    }
    else
    {
        rTimeControl.RestorePreviousTime();
        return rTimeControl.GetTimeStep() * 0.5;
    }
}
