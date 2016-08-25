#pragma once

namespace NuTo
{
//! @author Volker Hirthammer
//! @date July 09, 2015
//! @brief ...
class TimeDependencyBase
{
public:
    TimeDependencyBase();

    virtual ~TimeDependencyBase(){}

    virtual double GetTimeDependentFactor(double rTime) = 0;


};
} // namespace NuTo

