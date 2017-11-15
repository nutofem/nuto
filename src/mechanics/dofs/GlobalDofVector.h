#pragma once

#include "mechanics/dofs/DofVector.h"

namespace NuTo
{

class GlobalDofVector
{
public:
    NuTo::DofVector<double> J;
    NuTo::DofVector<double> K;

    GlobalDofVector& operator+=(const GlobalDofVector& rhs)
    {
        this->J += rhs.J;
        this->K += rhs.K;
        return *this;
    }
};

} /* NuTo */
