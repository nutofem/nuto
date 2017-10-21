//
// Created by phuschke on 6/19/17.
//

#pragma once

#include "mechanics/feti/FetiPreconditioner.h"

namespace NuTo
{
class FetiLumpedPreconditioner : public FetiPreconditioner
{
public:
    virtual void Compute(const StructureOutputBlockMatrix& hessian, const SparseMatrixType& B,
                         const std::map<int, int>&) override
    {
        const int numTotalDofs = B.cols();

        auto K = hessian.JJ.ExportToEigenSparseMatrix();
        K.conservativeResize(numTotalDofs, numTotalDofs);

        mLocalPreconditioner = B * K * B.transpose();
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
