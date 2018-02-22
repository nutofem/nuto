#pragma once

#include <vector>
#include "mechanics/timeIntegration/ExplicitRungeKutta.h"

namespace NuTo
{
namespace TimeIntegration
{
template <typename state>
class RK4 : public ExplicitRungeKutta<state>
{
public:
    RK4()
        : NuTo::TimeIntegration::ExplicitRungeKutta<state>({{}, {0.5}, {0.0, 0.5}, {0., 0., 1.}},
                                                           {1. / 6., 1. / 3., 1. / 3., 1. / 6.}, {0., 0.5, 0.5, 1.})

    {
    }
};
}
}
