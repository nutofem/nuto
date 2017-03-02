//
// Created by Thomas Titscher on 3/1/17.
//

#pragma once

#include "mechanics/dofSubMatrixSolvers/SolverBase.h"
#include "math/SparseDirectSolverPardiso.h"
#include "math/SparseMatrixCSR.h"

namespace NuTo
{

class SolverPardiso : public SolverBase
{
public:
    SolverPardiso(int rNumProcessors)
        : SolverBase()
        , mNumProcessors(rNumProcessors)
    {
    }
#ifdef HAVE_PARDISO
    virtual BlockFullVector<double> Solve(const BlockSparseMatrix& rMatrix,
                                          const BlockFullVector<double>& rVector) override
    {
        Eigen::VectorXd result;
        std::unique_ptr<NuTo::SparseMatrixCSR<double>> matrixForSolver = rMatrix.ExportToCSR();

        NuTo::SparseDirectSolverPardiso pardiso(mNumProcessors);
        pardiso.Solve(*matrixForSolver, rVector.Export(), result);

        return BlockFullVector<double>(result, rMatrix.GetDofStatus());
    }
#endif // HAVE_PARDISO
private:
    int mNumProcessors;
};
} // namespace NuTo
