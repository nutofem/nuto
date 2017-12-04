//
// Created by phuschke on 6/19/17.
//

#pragma once

#include <mpi.h>

#include "Eigen/Sparse"

#include "mechanics/structures/StructureOutputBlockMatrix.h"

namespace NuTo
{
class FetiPreconditioner
{
public:
    using SparseMatrixType = Eigen::SparseMatrix<double>;
    using VectorType = Eigen::VectorXd;

    virtual void Compute(const StructureOutputBlockMatrix& hessian, const SparseMatrixType& B,
                         const std::map<int, int>& lagrangeMultipliersGlobalIdToLocalId) = 0;

    virtual void AddScaling(const SparseMatrixType& scalingMatrix)
    {
        mLocalPreconditioner = scalingMatrix * mLocalPreconditioner * scalingMatrix;
    }

    //! \brief Applies the local preconditioner on the left and performs an MPI_Allreduce
    virtual VectorType ApplyOnTheLeft(const VectorType& x)
    {
        VectorType vec = mLocalPreconditioner * x;
        MPI_Allreduce(MPI_IN_PLACE, vec.data(), vec.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        return vec;
    }

protected:
    SparseMatrixType mLocalPreconditioner;
};
} // namespace NuTo
