//
// Created by phuschke on 6/19/17.
//

#pragma once

#include "mechanics/feti/FetiPreconditioner.h"
#include "mechanics/feti/ReverseMap.h"
#include <boost/range/adaptor/map.hpp>
template <class A, class B>
Eigen::SparseMatrix<double> ExtractSubMatrix(const Eigen::SparseMatrix<double>& mat, const A& rowIds, const B& colIds);

namespace NuTo
{
class FetiDirichletPreconditioner : public FetiPreconditioner
{
public:
    virtual void Compute(const StructureOutputBlockMatrix& hessian, const SparseMatrixType& B,
                         const std::map<int, int>& lagrangeMultipliersGlobalIdToLocalId) override
    {
        const int numTotalDofs = B.cols();

        std::set<int> internalDofIds;
        for (int j = 0; j < numTotalDofs; ++j)
        {
            internalDofIds.insert(j);
        }

        ReverseMap<int> reverseMap(lagrangeMultipliersGlobalIdToLocalId);
        std::vector<int> lagrangeMultiplierDofIds;
        for (auto& pair : reverseMap)
        {
            internalDofIds.erase(pair.first);
            lagrangeMultiplierDofIds.push_back(pair.first);
        }

        std::vector<int> internalDofIdsVec(internalDofIds.begin(), internalDofIds.end());

        SparseMatrixType H = hessian.ExportToEigenSparseMatrix();
        SparseMatrixType Kbb = ExtractSubMatrix(H, lagrangeMultiplierDofIds, lagrangeMultiplierDofIds);
        SparseMatrixType Kii = ExtractSubMatrix(H, internalDofIdsVec, internalDofIdsVec);
        SparseMatrixType Kbi = ExtractSubMatrix(H, lagrangeMultiplierDofIds, internalDofIdsVec);
        SparseMatrixType Kib = ExtractSubMatrix(H, internalDofIdsVec, lagrangeMultiplierDofIds);


        Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> luKii;
        luKii.compute(Kii);

        SparseMatrixType KiiInvTimesKib = luKii.solve(Kib);
        SparseMatrixType KbiTimesKiiInvTimesKib = Kbi * KiiInvTimesKib;
        SparseMatrixType Sbb = Kbb - KbiTimesKiiInvTimesKib;

        //
        //     | 0  0   |
        // S = | 0  Sbb |
        //
        SparseMatrixType S(numTotalDofs, numTotalDofs);


        for (size_t rowId = 0; rowId < lagrangeMultiplierDofIds.size(); ++rowId)
            for (size_t colId = 0; colId < lagrangeMultiplierDofIds.size(); ++colId)
                S.insert(lagrangeMultiplierDofIds[rowId], lagrangeMultiplierDofIds[colId]) = Sbb.coeff(rowId, colId);

        mLocalPreconditioner = B * S * B.transpose();
    }

    //! \brief Applies the local preconditioner on the left and performs an MPI_Allreduce
    VectorType ApplyOnTheLeft(const VectorType& x) override
    {
        assert(mLocalPreconditioner.cols() == x.rows());

        VectorType vec = mLocalPreconditioner * x;
        MPI_Allreduce(MPI_IN_PLACE, vec.data(), vec.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        return vec;
    }
};
} // namespace NuTo


template <class A, class B>
Eigen::SparseMatrix<double> ExtractSubMatrix(const Eigen::SparseMatrix<double>& mat, const A& rowIds, const B& colIds)
{
    Eigen::SparseMatrix<double> subMatrix(rowIds.size(), colIds.size());


    for (int k = 0; k < mat.outerSize(); ++k)
    {
        auto colIt = std::find(colIds.begin(), colIds.end(), k);
        if (colIt != colIds.end())
        {
            for (typename Eigen::SparseMatrix<double>::InnerIterator it(mat, k); it; ++it)
            {
                auto rowIt = std::find(rowIds.begin(), rowIds.end(), it.row());
                if (rowIt != rowIds.end())
                {
                    int rowId = std::distance(rowIds.begin(), rowIt);
                    int colId = std::distance(colIds.begin(), colIt);
                    subMatrix.insert(rowId, colId) = it.value();

                    ++rowId;
                }
            }
        }
    }

    subMatrix.makeCompressed();

    return subMatrix;
}
