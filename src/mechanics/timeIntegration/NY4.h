#pragma once

#include <vector>
#include "mechanics/timeIntegration/ExplicitNystroem.h"

namespace NuTo
{
namespace TimeIntegration
{
template <typename TState>
class NY4 : public ExplicitNystroem<TState>
{
public:
    NY4()
        : NuTo::TimeIntegration::ExplicitNystroem<TState>({{0., 0., 0.}, {1. / 8., 0., 0.}, {0., 0.5, 0.}},
                                                          {1. / 6, 2. / 3., 1. / 6.}, {1. / 6, 1. / 3., 0.},
                                                          {0., 1. / 2., 1.})
    {
    }
};
}
}
