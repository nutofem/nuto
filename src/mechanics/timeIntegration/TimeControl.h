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


    //! @brief Adjusts the timestep
    //! @param iterations: Number of iterations that were needed by the time integration scheme at the current time
    //! @param maxIterations: Maximum number of iterations allowed
    //! @param converged: Did the solution of the time integration scheme converge?
    void AdjustTimestep(int iterations, int maxIterations, bool converged);

    //! @brief Returnes if the time control has finished by reaching the final time
    //! @return true/false
    bool Finished (){return mCurrentTime >= mTimeFinal;}

    //! @brief Proceeds with the next time step    
    void Proceed();



    //! @brief Sets the timestep function that should be executed when the proceed function is called
    //! @param timestepFunction: function that should be executed when the proceed function is called
    void SetTimeStepFunction(std::function<double(TimeControl&, int, int, bool)> timestepFunction);

    //! @brief Resets the current time to the previous time
    void RestorePreviousTime(){mCurrentTime = mPreviousTime;}

    //! @brief Sets the timestep function to the default automatic timestepping method
    void UseDefaultAutomaticTimestepping();

    //! @brief Sets the timestep function to the default equidistant timestepping method
    void UseEquidistantTimestepping();

    //! @brief default automatic timestepping function that can be assigned to be the time stepping function of the time control
    //! @param timeControl: Reference to time control class
    //! @param iterations: Number of iterations that were needed by the time integration scheme at the current time
    //! @param maxIterations: Maximum number of iterations allowed
    //! @param converged: Did the solution of the time integration scheme converge?
    //! @return adjusted timestep
    static double DefaultAutomaticTimestepFunction(TimeControl& rTimeControl, int iterations, int maxIterations,bool converged);

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

    //! @brief Sets the time stepping to equidistant
    //! @param timestep: timestep value
    void SetTimeStep(double timeStep);

    //! @brief Sets the final time
    //! @param timefinal: final time
    void SetTimeFinal(double timeFinal);


    //temporary to remove all other time related members from timeIntegrationBase without bigger changes in derived classes solve routines
    // problem is that the postprocessor uses the timecontrol data, which is not fully implemented in all time integration schemes
    void SetCurrentTime(double curTime)
    {
        mCurrentTime = curTime;
    }


protected:



    double mCurrentTime             = 0.0;
    double mPreviousTime            = 0.0;
    double mTimeFinal               = 0.0;
    double mTimeStep                = 0.0;
    double mMinTimeStep             = 0.0;
    double mMaxTimeStep             = std::numeric_limits<double>::max();



#ifndef SWIG
    std::function<double(TimeControl&,int,int,bool)> mTimeStepFunction = [](TimeControl& timeControl, int iterations, int maxIterations,bool converged)->double{return timeControl.GetTimeStep();};
#endif
};

} // namespace NuTo
