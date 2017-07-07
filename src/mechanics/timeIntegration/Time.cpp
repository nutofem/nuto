#include "Time.h"



double NuTo::Time::Proceed()
{
// debug blocks ensure that the provided proceed funtion updates the time variables correctly
#ifndef NDEBUG
    double currentTimeBeforeProceed = mCurrentTime;
#endif

    mCurrentTime = mTimestepFunction(mCurrentTime);

#ifndef NDEBUG
    if(currentTimeBeforeProceed == mCurrentTime)
        throw MechanicsException(__PRETTY_FUNCTION__,"The current proceed function of the class does not increase the current time!");
#endif

    return mCurrentTime;
}

void NuTo::Time::SetEquidistantTimestepping(double timestep)
{
    if (timestep<=0)
        throw MechanicsException(__PRETTY_FUNCTION__,"Timestep must be a positive number!");

    mTimestepFunction =  [timestep](double curTime)->double
                        {
                            return curTime + timestep;
                        };
}

void NuTo::Time::SetTimestepFunction(std::function<double (double)> timestepFunction)
{
    mTimestepFunction = timestepFunction;
}
