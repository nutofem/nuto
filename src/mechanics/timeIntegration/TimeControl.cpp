#include "TimeControl.h"



void NuTo::TimeControl::ScaleTimeStep(double scaleFactor)
{
    if (scaleFactor<= 0.0)
        throw MechanicsException(__PRETTY_FUNCTION__,"Scaling factor must be a positive number");

    mTimeStepScaleFactor *= scaleFactor;
    UpdateTimeStep();
}

double NuTo::TimeControl::Proceed()
{
    mPreviousTime = mCurrentTime;

    UpdateTimeStep();
    mCurrentTime += mTimeStep;

    return mCurrentTime;
}



void NuTo::TimeControl::SetTimeStep(double TimeStep)
{
    if (TimeStep<=0)
        throw MechanicsException(__PRETTY_FUNCTION__,"Timestep must be a positive number!");

    SetTimeStepFunction([TimeStep]()->double
                        {
                            return TimeStep;
                        });
}



void NuTo::TimeControl::SetTimeStepFunction(std::function<double ()> TimeStepFunction)
{
    mTimeStepFunction       = TimeStepFunction;
    ResetTimeStepScaleFactor();
    UpdateTimeStep();
}

void NuTo::TimeControl::SetMaxTimeStep(double rMaxTimeStep)
{
    if (rMaxTimeStep <=0.0)
        throw MechanicsException(__PRETTY_FUNCTION__,"Maximal timestep must be a positive number!");
    if (rMaxTimeStep<mMinTimeStep)
        throw MechanicsException(__PRETTY_FUNCTION__,"Maximal timestep must be bigger than minimal Timestep!");
    mMaxTimeStep = rMaxTimeStep;
}

void NuTo::TimeControl::SetMinTimeStep(double rMinTimeStep)
{
    if (rMinTimeStep>mMaxTimeStep)
        throw MechanicsException(__PRETTY_FUNCTION__,"Minimal timestep must be smaller than maximal Timestep!");
    mMinTimeStep = rMinTimeStep;
}

void NuTo::TimeControl::UpdateTimeStep()
{
    mTimeStep = mTimeStepScaleFactor * mTimeStepFunction();

    if(mTimeStep > mMaxTimeStep)
        mTimeStep =  mMaxTimeStep;

    if (mTimeStep <= 0.0)
        throw MechanicsException(__PRETTY_FUNCTION__,"Current timestep is 0 or negative!");
    if (mTimeStep<mMinTimeStep)
        throw MechanicsException(__PRETTY_FUNCTION__,"Current timestep is lower than minimum Timestep!");
}


