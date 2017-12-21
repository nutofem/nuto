#pragma once

#include "mechanics/dofs/DofCalcContainer.h"

namespace NuTo
{
template <typename T>
using DofVector = DofCalcContainer<Eigen::Matrix<T, Eigen::Dynamic, 1>>;
} /* NuTo */
