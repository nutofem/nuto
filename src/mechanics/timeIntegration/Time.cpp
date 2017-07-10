#include "Time.h"

using namespace NuTo;

void Time::ScaleTimestep(double scaleFactor)
{
    if (scaleFactor <= 0.0)
        throw MechanicsException(__PRETTY_FUNCTION__, "Scaling factor must be a positive number");
    mTimestepScaleFactor *= scaleFactor;

    if (GetTimestep() < mMinTimeStep)
        throw MechanicsException(__PRETTY_FUNCTION__, "Current timestep is lower than Minimum!");
}


double Time::Proceed()
{
    mPreviousTime = mCurrentTime;
    mCurrentTime += GetTimestep();
    if (mPreviousTime == mCurrentTime)
    {
        if (mTimestepScaleFactor > 0.0)
            throw MechanicsException(__PRETTY_FUNCTION__,
                                     "The current proceed function of the class does not increase the current time!");
        else
            throw MechanicsException(__PRETTY_FUNCTION__, "Timestep scaling factor 0 or negative!");
    }
    return mCurrentTime;
}


void Time::SetEquidistantTimestepping(double timestep)
{
    if (timestep <= 0)
        throw MechanicsException(__PRETTY_FUNCTION__, "Timestep must be a positive number!");

    mTimestepFunction = [timestep](double curTime) -> double { return timestep; };
}


void Time::SetTimestepFunction(std::function<double(double)> timestepFunction)
{
    mTimestepFunction = timestepFunction;
}


double Time::GetTimestep() const
{
    double currentTimestep = mTimestepScaleFactor * mTimestepFunction(mCurrentTime);
    if (currentTimestep > mMaxTimeStep)
        return mMaxTimeStep;

    return currentTimestep;
}
