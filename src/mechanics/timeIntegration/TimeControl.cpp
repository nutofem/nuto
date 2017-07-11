#include "TimeControl.h"



void NuTo::TimeControl::ScaleTimestep(double scaleFactor)
{
    if (scaleFactor<= 0.0)
        throw MechanicsException(__PRETTY_FUNCTION__,"Scaling factor must be a positive number");

    mTimestepScaleFactor *= scaleFactor;
    UpdateTimestep();
}

double NuTo::TimeControl::Proceed()
{
    mPreviousTime = mCurrentTime;

    UpdateTimestep();
    mCurrentTime += mTimestep;

    return mCurrentTime;
}



void NuTo::TimeControl::SetEquidistantTimestepping(double timestep)
{
    if (timestep<=0)
        throw MechanicsException(__PRETTY_FUNCTION__,"Timestep must be a positive number!");

    SetTimestepFunction([timestep]()->double
                        {
                            return timestep;
                        });
}



void NuTo::TimeControl::SetTimestepFunction(std::function<double ()> timestepFunction)
{
    mTimestepFunction       = timestepFunction;
    ResetTimestepScaleFactor();
    UpdateTimestep();
}

void NuTo::TimeControl::SetMaxTimestep(double rMaxTimeStep)
{
    if (rMaxTimeStep <=0.0)
        throw MechanicsException(__PRETTY_FUNCTION__,"Maximal timestep must be a positive number!");
    if (rMaxTimeStep<mMinTimeStep)
        throw MechanicsException(__PRETTY_FUNCTION__,"Maximal timestep must be bigger than minimal timestep!");
    mMaxTimeStep = rMaxTimeStep;
}

void NuTo::TimeControl::SetMinTimestep(double rMinTimeStep)
{
    if (rMinTimeStep>mMaxTimeStep)
        throw MechanicsException(__PRETTY_FUNCTION__,"Minimal timestep must be smaller than maximal timestep!");
    mMinTimeStep = rMinTimeStep;
}

void NuTo::TimeControl::UpdateTimestep()
{
    mTimestep = mTimestepScaleFactor * mTimestepFunction();

    if(mTimestep > mMaxTimeStep)
        mTimestep =  mMaxTimeStep;

    if (mTimestep <= 0.0)
        throw MechanicsException(__PRETTY_FUNCTION__,"Current timestep is 0 or negative!");
    if (mTimestep<mMinTimeStep)
        throw MechanicsException(__PRETTY_FUNCTION__,"Current timestep is lower than minimum timestep!");
}


