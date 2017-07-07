#pragma once

#include "mechanics/MechanicsException.h"

#include <functional>

namespace NuTo
{


class Time
{
public:

    //! @brief ctor
    Time() = default;

    //! @brief copy ctor
    Time(const Time& rOther) = default;

    //! @brief move ctor
    Time(Time&& rOther) = default;

    //! @brief copy assignment
    Time& operator=(const Time& rOther) = default;

    //! @brief move assignment
    Time& operator=(Time&& rOther) = default;

    //! @brief destructor
    ~Time() = default;



    //! @brief Proceeds with the next time step
    //! @return current time value
    double Proceed();

    //! @brief Sets the time stepping to equidistant
    //! @param timestep: timestep value
    void SetEquidistantTimestepping(double timestep);

    //! @brief Sets the timestep function that should be executed when the proceed function is called
    //! @param timestepFunction: function that should be executed when the proceed function is called
    void SetTimestepFunction(std::function<double(double)> timestepFunction);


    // Getter
    // ------
    //! @brief Gets the current time
    //! @return current time
    double GetCurrentTime() const {return mCurrentTime;}


protected:

    double mCurrentTime     = 0;

    std::function<double(double)> mTimestepFunction = [](double curTime)->double{throw MechanicsException(__PRETTY_FUNCTION__,"No timestepping method selected!");};
};

} // namespace NuTo
