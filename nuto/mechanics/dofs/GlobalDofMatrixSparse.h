#pragma once

#include "nuto/mechanics/dofs/DofMatrixSparse.h"

namespace NuTo
{
class GlobalDofMatrixSparse
{
public:
    DofMatrixSparse<double> JJ;
    DofMatrixSparse<double> JK;
    DofMatrixSparse<double> KJ;
    DofMatrixSparse<double> KK;

    GlobalDofMatrixSparse& operator+=(const GlobalDofMatrixSparse& rhs)
    {
        this->JJ += rhs.JJ;
        this->JK += rhs.JK;
        this->KJ += rhs.KJ;
        this->KK += rhs.KK;
        return *this;
    }
};
} /* NuTo */
