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
    SolverPardiso(int rNumProcessors, bool rShowTime = true)
        : SolverBase()
#ifdef HAVE_PARDISO
        , mNumProcessors(rNumProcessors)
        , mShowTime(rShowTime)
#endif // HAVE_PARDISO
    {
    }
#ifdef HAVE_PARDISO
    virtual BlockFullVector<double> Solve(const BlockSparseMatrix& rMatrix,
                                          const BlockFullVector<double>& rVector) override
    {
        Eigen::VectorXd result;
        std::unique_ptr<NuTo::SparseMatrixCSR<double>> matrixForSolver = rMatrix.ExportToCSR();
        matrixForSolver->SetOneBasedIndexing();

        int verboseLevel = mShowTime ? 1 : 0;
        NuTo::SparseDirectSolverPardiso pardiso(mNumProcessors, verboseLevel);
        pardiso.SetShowTime(mShowTime);
        pardiso.Solve(*matrixForSolver, rVector.Export(), result);

        return BlockFullVector<double>(result, rMatrix.GetDofStatus());
    }

private:
    int mNumProcessors;
    bool mShowTime;
#endif // HAVE_PARDISO
};
} // namespace NuTo
