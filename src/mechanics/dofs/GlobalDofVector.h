#pragma once

#include "mechanics/dofs/DofVector.h"

namespace NuTo
{

class GlobalDofVector
{
public:
    NuTo::DofVector<double> J;
    NuTo::DofVector<double> K;
};

} /* NuTo */
