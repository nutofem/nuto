#pragma once

#include "mechanics/dofs/DofMatrixSparse.h"

namespace NuTo
{
class GlobalDofMatrixSparse
{
public:
    DofMatrixSparse<double> JJ;
    DofMatrixSparse<double> JK;
    DofMatrixSparse<double> KJ;
    DofMatrixSparse<double> KK;
};
} /* NuTo */
