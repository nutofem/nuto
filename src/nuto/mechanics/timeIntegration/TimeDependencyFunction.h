#ifndef TIMEDEPENDENCYFUNCTION
#define TIMEDEPENDENCYFUNCTION

#include <functional>

namespace NuTo
{
//! @author Volker Hirthammer
//! @date July 09, 2015
//! @brief ...
class TimeDependencyFunction : public TimeDependencyBase
{
public:
    TimeDependencyFunction(const std::function<double (double rTime)>& rTimeDependencyFunction)
        : TimeDependencyBase(),
          mTimeDependencyFunction(rTimeDependencyFunction)
    {}

    virtual ~TimeDependencyFunction(){}

    virtual double GetTimeDependentFactor(double rTime) override
    {
        return mTimeDependencyFunction(rTime);
    }


private:
    std::function<double (double rTime)> mTimeDependencyFunction;
};
} // namespace NuTo


#endif // TIMEDEPENDENCYFUNCTION

