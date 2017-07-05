//
// Created by phuschke on 3/1/17.
//

#pragma once

#include "mechanics/dofSubMatrixSolvers/SolverBase.h"

namespace NuTo
{
template <class Solver>
class SolverEigen : public SolverBase
{
public:
    SolverEigen()
        : SolverBase()
    {
    }
    virtual BlockFullVector<double> Solve(const BlockSparseMatrix& rMatrix,
                                          const BlockFullVector<double>& rVector) override
    {
        Eigen::VectorXd result;
        Solver solver;
        solver.compute(rMatrix.ExportToEigenSparseMatrix());
        return BlockFullVector<double>(solver.solve(rVector.Export()), rMatrix.GetDofStatus());
    }
};
} // namespace NuTo
