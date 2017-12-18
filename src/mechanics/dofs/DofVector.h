#pragma once

#include "mechanics/dofs/DofCalcContainer.h"

namespace NuTo
{
template <typename T>
using DofVector = DofCalcContainer<Eigen::Matrix<T, Eigen::Dynamic, 1>>;

//! @brief export the dofs entries of a DofVector to a Eigen::VectorXd
//! @tparam T numeric type
//! @param v dof vector to export
//! @param dofs dof types to export
//! @return continuous vector containing the combined subvectors of v
template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> ToEigen(const DofVector<T>& v, std::vector<DofType> dofs)
{
    int totalRows = 0;
    for (auto dof : dofs)
        totalRows += v[dof].rows();
    Eigen::Matrix<T, Eigen::Dynamic, 1> combined(totalRows);

    int currentStartRow = 0;
    for (auto dof : dofs)
    {
        const int rows = v[dof].rows();
        combined.segment(currentStartRow, rows) = v[dof];
        currentStartRow += rows;
    }
    return combined;
}

} /* NuTo */
