#include "TimeControl.h"



void NuTo::TimeControl::ScaleTimeStep(double scaleFactor)
{
    if (scaleFactor<= 0.0)
        throw MechanicsException(__PRETTY_FUNCTION__,"Scaling factor must be a positive number");

    mTimeStepScaleFactor *= scaleFactor;
    UpdateTimeStep(-1.,-1.,true);
}

void NuTo::TimeControl::Proceed()
{
    Proceed(-1.,-1.,true);
}

void NuTo::TimeControl::Proceed(double iterations, double maxIterations, bool convergence)
{
    if(mCurrentTime == mTimeFinal)
    {
        mFinished = true;
        return;
    }

    UpdateTimeStep(iterations,maxIterations,convergence);
    if(mCurrentTime < mTimeFinal)
    {
        if(!convergence && mPreviousTime !=mCurrentTime)
            throw MechanicsException(__PRETTY_FUNCTION__, "No convergence with the current maximum number of "
                                                          "iterations, either use automatic time stepping, "
                                                          "reduce the time step or the minimal line search cut "
                                                          "back factor. In case you provided a custom timestepping function "
                                                          "and intend to use some kind of automatic timestepping, call "
                                                          "the RestorePreviosTime() function of the time control before reducing "
                                                          "the timestep.");
        mPreviousTime = mCurrentTime;
        mCurrentTime += mTimeStep;
        if (mCurrentTime > mTimeFinal)
            mCurrentTime = mTimeFinal;
    }
}



void NuTo::TimeControl::SetTimeStep(double TimeStep)
{
    if (TimeStep<=0)
        throw MechanicsException(__PRETTY_FUNCTION__,"Timestep must be a positive number!");

    mTimeStep = TimeStep;
}

void NuTo::TimeControl::SetTimeFinal(double timefinal)
{
    if (timefinal <= mCurrentTime)
    {
        throw MechanicsException(__PRETTY_FUNCTION__,"Final time must be larger than current time!");
    }
    mTimeFinal          = timefinal;
    mFinished           = false;
}



void NuTo::TimeControl::SetTimeStepFunction(std::function<double(double, double, bool)> TimeStepFunction)
{
    mTimeStepFunction       = TimeStepFunction;
}

void NuTo::TimeControl::UseDefaultAutomaticTimestepping()
{
    SetTimeStepFunction([this](double iterations, double maxIterations, bool convergence)->double
                        {
                            if(iterations < 0.25 * maxIterations)
                                return mTimeStep * 1.5;
                            if(!convergence)
                                RestorePreviosTime();
                                return mTimeStep * 0.5;
                            return mTimeStep;
                        });
}

void NuTo::TimeControl::UseEquidistantTimestepping()
{
    SetTimeStepFunction([this](double iterations, double maxIterations, bool convergence)->double
                        {
                            return mTimeStep;
                        });
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

void NuTo::TimeControl::UpdateTimeStep(double iterations, double maxIterations, bool convergence)
{
//    mTimeStep = mTimeStepScaleFactor * mTimeStepFunction(iterations,maxIterations);
    mTimeStep = mTimeStepFunction(iterations,maxIterations,convergence);

    if(mTimeStep > mMaxTimeStep)
        mTimeStep =  mMaxTimeStep;

    if (mTimeStep <= 0.0)
        throw MechanicsException(__PRETTY_FUNCTION__,"Current timestep is 0 or negative!");
    if (mTimeStep<mMinTimeStep)
        throw MechanicsException(__PRETTY_FUNCTION__,"Current timestep is lower than minimum Timestep!");
}


