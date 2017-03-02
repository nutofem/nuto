//
// Created by phuschke on 3/1/17.
//


#pragma once

#include "mechanics/dofSubMatrixSolvers/SolverBase.h"
#include "math/SparseDirectSolverMUMPS.h"
#include "math/SparseMatrixCSR.h"

namespace NuTo
{

    class SolverMUMPS : public SolverBase
    {
    public:
        SolverMUMPS() : SolverBase()
        {}
        virtual BlockFullVector<double> Solve(const BlockSparseMatrix& rMatrix, const BlockFullVector<double>& rVector) override
        {

            Eigen::VectorXd result;
            std::unique_ptr<NuTo::SparseMatrixCSR<double>> matrixForSolver = rMatrix.ExportToCSR();
            matrixForSolver->SetOneBasedIndexing();

            NuTo::SparseDirectSolverMUMPS mySolver;

            mySolver.Solve(*matrixForSolver, rVector.Export(), result);

            return BlockFullVector<double>(result, rMatrix.GetDofStatus());

        }
    };
} // namespace NuTo
