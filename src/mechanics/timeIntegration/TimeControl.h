#pragma once

#include "mechanics/MechanicsException.h"

#include <functional>
#include <limits>

namespace NuTo
{


class TimeControl
{
public:

    //! @brief ctor
    TimeControl() = default;

    //! @brief copy ctor
    TimeControl(const TimeControl& rOther) = default;

    //! @brief move ctor
    TimeControl(TimeControl&& rOther) = default;

    //! @brief copy assignment
    TimeControl& operator=(const TimeControl& rOther) = delete;

    //! @brief move assignment
    TimeControl& operator=(TimeControl&& rOther) = delete;

    //! @brief destructor
    ~TimeControl() = default;


    //! @brief Returnes if the time control has finiched by reaching the final time
    //! @return true/false
    bool Finished (){return mFinished;}

    //! @brief Scales the timestep by the provided factor
    //! @param scaleFactor: scaling factor (<1 decrease and >1 increase)
    void ScaleTimeStep(double scaleFactor);

    //! @brief Proceeds with the next time step
    //! @return current time value
    void Proceed();

    //! @brief Proceeds with the next time step
    //! @return current time value
    void Proceed(double iterations, double maxIterations, bool convergence);

    //! @brief Sets the time stepping to equidistant
    //! @param timestep: timestep value
    void SetTimeStep(double timestep);

    //! @brief Sets the final time
    //! @param timefinal: final time
    void SetTimeFinal(double timefinal);

    //! @brief Sets the timestep function that should be executed when the proceed function is called
    //! @param timestepFunction: function that should be executed when the proceed function is called
    void SetTimeStepFunction(std::function<double(double, double,bool)> timestepFunction);

    //! @brief Resets the timestep scaling factor to 1 wich means the timestep is exactly as provided by the timestepfunction
    void ResetTimeStepScaleFactor()
    {
        mTimeStepScaleFactor    = 1.0;
    }

    //! @brief Resets the current time to the previous time
    void RestorePreviosTime()
    {
        mCurrentTime = mPreviousTime;
    }

    //! @brief Sets the timestep function to the default automatic timestepping method
    void UseDefaultAutomaticTimestepping();

    //! @brief Sets the timestep function to the default equidistant timestepping method
    void UseEquidistantTimestepping();

    // Getter
    // ------
    //! @brief Gets the current time
    //! @return current time
    double GetCurrentTime() const {return mCurrentTime;}

    //! @brief Gets the previous time
    //! @return previous time
    double GetPreviousTime() const {return mPreviousTime;}

    //! @brief Gets the timestep
    //! @return timestep
    double GetTimeStep() const {return mTimeStep;}

    //! @brief Gets the minimal timestep
    //! @return minimal timestep
    double GetMinTimeStep() const {return mMinTimeStep;}

    //! @brief Gets the maximum timestep
    //! @return maximum timestep
    double GetMaxTimeStep() const {return mMaxTimeStep;}

    // Setter
    // ------

    //! @brief sets the maximum time step for the time integration procedure
    void SetMaxTimeStep(double rMaxTimeStep);

    //! @brief sets the minimum time step for the time integration procedure
    void SetMinTimeStep(double rMinTimeStep);


    //temporary to remove all other time related members from timeIntegrationBase without bigger changes in derived classes solve routines
    void SetCurrentTime(double curTime)
    {
        mCurrentTime = curTime;
    }

protected:

    void UpdateTimeStep(double iterations, double maxIterations, bool convergence);



    double mCurrentTime             = 0.0;
    double mPreviousTime            = 0.0;
    double mTimeFinal               = 0.0;
    double mTimeStepScaleFactor     = 1.0;
    double mTimeStep                = 0.0;
    double mMinTimeStep             = 0.0;
    double mMaxTimeStep             = std::numeric_limits<double>::max();

//    bool mAutomaticTimestepping     = false;
    bool mFinished                  = false;


#ifndef SWIG
    std::function<double(double,double,bool)> mTimeStepFunction = [this](double iterations, double maxIterations,bool convergence)->double{return mTimeStep;};
#endif
};

} // namespace NuTo
