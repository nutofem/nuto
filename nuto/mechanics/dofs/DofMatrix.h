#pragma once

#include "nuto/mechanics/dofs/DofMatrixContainer.h"

namespace NuTo
{
//! @brief dof container that is also capable of performing calculations.
template <typename T>
using DofMatrix = DofMatrixContainer<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>;
}

