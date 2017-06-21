//
// Created by phuschke on 6/19/17.
//

#pragma once

#include "mechanics/feti/FetiPreconditioner.h"

namespace NuTo
{
class FetiIdentityPreconditioner : public FetiPreconditioner
{
public:

    virtual void Compute(const StructureOutputBlockMatrix& hessian, const SparseMatrixType& B, const std::vector<int> lagrangeMultiplierDofIds) override
    {}

    //! \brief Applies the local preconditioner on the left and performs an MPI_Allreduce
    VectorType ApplyOnTheLeft(const VectorType& x) override
    {
        return x;
    }
};
} // namespace NuTo
