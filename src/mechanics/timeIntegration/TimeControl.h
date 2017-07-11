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


    //! @brief Scales the timestep by the provided factor
    //! @param scaleFactor: scaling factor (<1 decrease and >1 increase)
    void ScaleTimestep(double scaleFactor);

    //! @brief Proceeds with the next time step
    //! @return current time value
    double Proceed();

    //! @brief Sets the time stepping to equidistant
    //! @param timestep: timestep value
    void SetEquidistantTimestepping(double timestep);

    //! @brief Sets the timestep function that should be executed when the proceed function is called
    //! @param timestepFunction: function that should be executed when the proceed function is called
    void SetTimestepFunction(std::function<double()> timestepFunction);

    //! @brief Resets the timestep scaling factor to 1 wich means the timestep is exactly as provided by the timestepfunction
    void ResetTimestepScaleFactor()
    {
        mTimestepScaleFactor    = 1.0;
    }

    //! @brief Resets the current time to the previous time
    void RestorePreviosTime()
    {
        mCurrentTime = mPreviousTime;
    }

    // Getter
    // ------
    //! @brief Gets the current time
    //! @return current time
    double GetCurrentTime() const {return mCurrentTime;}

    //! @brief Gets the timestep
    //! @return timestep
    double GetTimestep() const
    {
        return mTimestep;
    }

    // Setter
    // ------

    //! @brief sets the maximum time step for the time integration procedure
    void SetMaxTimestep(double rMaxTimeStep);

    //! @brief sets the minimum time step for the time integration procedure
    void SetMinTimestep(double rMinTimeStep);


    //temporary to remove all other time related members from timeIntegrationBase without bigger changes in derived classes solve routines
    void SetCurrentTime(double curTime)
    {
        mCurrentTime = curTime;
    }

protected:

    void UpdateTimestep();



    double mCurrentTime             = 0.0;
    double mPreviousTime            = 0.0;
    double mTimestepScaleFactor     = 1.0;
    double mTimestep                = 0.0;
    double mMinTimeStep             = 0.0;
    double mMaxTimeStep             = std::numeric_limits<double>::max();

    std::function<double()> mTimestepFunction = []()->double{throw MechanicsException(__PRETTY_FUNCTION__,"No timestepping method selected!");};
};

} // namespace NuTo
