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
    SolverMUMPS(bool rShowTime = true)
        : SolverBase()
        , mShowTime(rShowTime)
    {
    }
    virtual BlockFullVector<double> Solve(const BlockSparseMatrix& rMatrix,
                                          const BlockFullVector<double>& rVector) override
    {

        Eigen::VectorXd result;
        std::unique_ptr<NuTo::SparseMatrixCSR<double>> matrixForSolver = rMatrix.ExportToCSR();
        matrixForSolver->SetOneBasedIndexing();

        NuTo::SparseDirectSolverMUMPS solver;
        solver.SetShowTime(mShowTime);

        solver.Solve(*matrixForSolver, rVector.Export(), result);

        return BlockFullVector<double>(result, rMatrix.GetDofStatus());
    }

private:
    bool mShowTime;
};
} // namespace NuTo
