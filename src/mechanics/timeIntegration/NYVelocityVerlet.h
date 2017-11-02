#pragma once

#include <vector>
#include "mechanics/timeIntegration/ExplicitNystroem.h"

namespace NuTo
{
namespace TimeIntegration
{
template <typename TState>
class NYVelocityVerlet : public ExplicitNystroem<TState>
{
public:
    NYVelocityVerlet()
        : NuTo::TimeIntegration::ExplicitNystroem<TState>({{0., 0.}, {0.5, 0.}}, {0.5, 0.5}, {0.5, 0.},
                                                                    {0., 1.})
    {
    }
};
}
}
