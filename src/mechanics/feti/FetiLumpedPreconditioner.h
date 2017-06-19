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
    //! \brief Applies the local preconditioner on the left and performs an MPI_Allreduce
    VectorType ApplyOnTheLeft(const VectorType& x) override
    {
        return x;
    }
};
} // namespace NuTo
