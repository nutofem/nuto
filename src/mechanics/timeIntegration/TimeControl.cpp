#include "TimeControl.h"



void NuTo::TimeControl::Proceed()
{
    AdjustTimestep(0.,0.,true);
}

void NuTo::TimeControl::AdjustTimestep(int iterations, int maxIterations, bool converged)
{
    if(mCurrentTime == mTimeFinal)
    {
        mFinished = true;
        return;
    }

    UpdateTimeStep(iterations,maxIterations,converged);

    if(!converged && mPreviousTime !=mCurrentTime)
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



void NuTo::TimeControl::SetTimeStepFunction(std::function<double(TimeControl&,int, int, bool)> TimeStepFunction)
{
    mTimeStepFunction       = TimeStepFunction;
}

void NuTo::TimeControl::UseDefaultAutomaticTimestepping()
{
    SetTimeStepFunction([this](TimeControl& timeControl, int iterations, int maxIterations, bool converged)->double
                        {
                            if(iterations < 0.25 * maxIterations)
                                return timeControl.GetTimeStep() * 1.5;
                            if(!converged)
                                RestorePreviosTime();
                                return timeControl.GetTimeStep() * 0.5;
                            return timeControl.GetTimeStep();
                        });
}

void NuTo::TimeControl::UseEquidistantTimestepping()
{
    SetTimeStepFunction([this](TimeControl& timeControl, int iterations, int maxIterations, bool converged)->double
                        {
                            return timeControl.GetTimeStep();
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

void NuTo::TimeControl::UpdateTimeStep(int iterations, int maxIterations, bool converged)
{
    mTimeStep = mTimeStepFunction(*this,iterations,maxIterations,converged);

    if(mTimeStep > mMaxTimeStep)
        mTimeStep =  mMaxTimeStep;

    if (mTimeStep <= 0.0)
        throw MechanicsException(__PRETTY_FUNCTION__,"Current timestep is 0 or negative!");
    if (mTimeStep<mMinTimeStep)
        throw MechanicsException(__PRETTY_FUNCTION__,"Current timestep is lower than minimum Timestep!");
}


