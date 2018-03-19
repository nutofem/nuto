#pragma once

#include "nuto/mechanics/dofs/DofVector.h"

namespace NuTo
{
template <typename T>
inline int TotalRows(const DofVector<T>& v, const std::vector<DofType>& dofs)
{
    int totalRows = 0;
    for (auto dof : dofs)
        totalRows += v[dof].rows();
    return totalRows;
}

//! @brief export the dofs entries of a DofVector to a Eigen::VectorXT
//! @tparam T numeric type
//! @param v dof vector to export
//! @param dofs dof types to export
//! @return continuous vector containing the combined subvectors of v
template <typename T>
inline Eigen::Matrix<T, Eigen::Dynamic, 1> ToEigen(const DofVector<T>& v, std::vector<DofType> dofs)
{
    Eigen::Matrix<T, Eigen::Dynamic, 1> combined(TotalRows(v, dofs));

    int currentStartRow = 0;
    for (auto dof : dofs)
    {
        const int rows = v[dof].rows();
        combined.segment(currentStartRow, rows) = v[dof];
        currentStartRow += rows;
    }
    return combined;
}

//! @brief imports a values into a properly sized DofVector
//! @param source eigen vector whose values are imported
//! @param dofs dof types to import
//! @param rDestination properly sized dof vector
template <typename T>
inline void FromEigen(const Eigen::Matrix<T, Eigen::Dynamic, 1>& source, std::vector<DofType> dofs,
                      DofVector<T>* rDestination)
{
    assert(rDestination != nullptr);
    assert(TotalRows(*rDestination, dofs) == source.rows());
    int currentStartRow = 0;
    for (auto dof : dofs)
    {
        const int rows = (*rDestination)[dof].rows();
        (*rDestination)[dof] = source.segment(currentStartRow, rows);
        currentStartRow += rows;
    }
}
} /* NuTo */
