#pragma once

#include "mechanics/dofs/DofMatrixSparse.h"

namespace NuTo
{
template <typename T>
inline int TotalRows(const DofMatrixSparse<T>& v, std::vector<DofType> dofs)
{
    int total = 0;
    for (auto dof : dofs)
        total += v(dof, dofs.front()).rows();
    return total;
}

template <typename T>
inline int TotalCols(const DofMatrixSparse<T>& v, std::vector<DofType> dofs)
{
    int total = 0;
    for (auto dof : dofs)
        total += v(dofs.front(), dof).cols();
    return total;
}

template <typename T>
inline int TotalNonZeros(const DofMatrixSparse<T>& v, std::vector<DofType> dofs)
{
    int total = 0;
    for (auto dof0 : dofs)
        for (auto dof1 : dofs)
            total += v(dof0, dof1).nonZeros();
    return total;
}

//! @brief export the dofs entries of a DofMatrixSparse to a Eigen::SparseMatrix
//! @tparam T numeric type
//! @param v dof matrix to export
//! @param dofs dof types to export
//! @return continuous matrix containing the combined submatrices of v
template <typename T>
inline Eigen::SparseMatrix<T> ToEigen(const DofMatrixSparse<T>& v, std::vector<DofType> dofs)
{
    std::vector<Eigen::Triplet<T>> triplets;
    int startRow = 0;
    for (auto dofRow : dofs)
    {
        int startCol = 0;
        for (auto dofCol : dofs)
        {
            const auto& mat = v(dofRow, dofCol);
            for (int k = 0; k < mat.outerSize(); ++k)
                for (Eigen::SparseMatrix<double>::InnerIterator it(mat, k); it; ++it)
                    triplets.emplace_back(Eigen::Triplet<T>(it.row() + startRow, it.col() + startCol, it.value()));

            startCol += mat.cols();
        }
        startRow += v(dofRow, dofRow).rows();
    }
    Eigen::SparseMatrix<T> combined(TotalRows(v, dofs), TotalCols(v, dofs));
    combined.setFromTriplets(triplets.begin(), triplets.end());
    combined.makeCompressed();
    return combined;
}
}
