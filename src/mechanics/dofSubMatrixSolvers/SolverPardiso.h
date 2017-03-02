//
// Created by phuschke on 3/2/17.
//

//
// Created by phuschke on 3/1/17.
//


#pragma once

#include "mechanics/dofSubMatrixSolvers/SolverBase.h"
#include "math/SparseDirectSolverPardiso.h"
#include "mechanics/dofSubMatrixStorage/BlockSparseMatrix.h"
#include "mechanics/dofSubMatrixStorage/BlockFullVector.h"
#include "math/SparseMatrixCSR.h"
#include <eigen3/Eigen/Dense>

namespace NuTo
{

    class SolverPardiso : public SolverBase
    {
    public:
        SolverPardiso() : SolverBase()
        {}

        virtual BlockFullVector<double> Solve(const BlockSparseMatrix& rMatrix, const BlockFullVector<double>& rVector) override
        {
            throw MathException(__PRETTY_FUNCTION__, "not implemented");
        }
    };
} // namespace NuTo
