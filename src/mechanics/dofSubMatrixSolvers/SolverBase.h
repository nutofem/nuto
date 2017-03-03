//
// Created by phuschke on 3/1/17.
//

#pragma once

#include "mechanics/dofSubMatrixStorage/BlockSparseMatrix.h"
#include "mechanics/dofSubMatrixStorage/BlockFullVector.h"
#include <eigen3/Eigen/Dense>

namespace NuTo
{
    class SolverBase
    {
    public:
        virtual ~SolverBase() = default;

        virtual BlockFullVector<double> Solve(const BlockSparseMatrix& rMatrix, const BlockFullVector<double>& rVector) = 0;
    };
} // namespace NuTo