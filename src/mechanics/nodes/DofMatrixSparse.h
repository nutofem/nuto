#pragma once

#include "mechanics/nodes/DofMatrixContainer.h"
#include <eigen3/Eigen/Sparse>

namespace NuTo
{
//! @brief dof container that is also capable of performing calculations.
template <typename T>
using DofMatrixSparse = DofMatrixContainer<Eigen::SparseMatrix<T>>;
}

