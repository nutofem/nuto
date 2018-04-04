#pragma once
#include <Eigen/Core>
#include "nuto/mechanics/dofs/DofVector.h"
#include "nuto/mechanics/dofs/DofMatrixSparse.h"
#include "nuto/mechanics/constraints/Constraints.h"

namespace NuTo
{

DofVector<double> Solve(const DofMatrixSparse<double>& K, const DofVector<double>& f, Constraint::Constraints& bcs,
                        std::vector<DofType> dofs, std::string solver = "EigenSparseLU");

DofVector<double> SolveTrialState(const DofMatrixSparse<double>& K, const DofVector<double>& f, double oldTime,
                                  double newTime, Constraint::Constraints& bcs, std::vector<DofType> dofs,
                                  std::string solver = "EigenSparseLU");

class ConstrainedSystemSolver
{
public:
    ConstrainedSystemSolver(Constraint::Constraints& bcs, std::vector<DofType> dofs, std::string solver);
    DofVector<double> Solve(const DofMatrixSparse<double>& A, const DofVector<double>& b);
    DofVector<double> SolveTrialState(const DofMatrixSparse<double>& A, const DofVector<double>& b, double oldTime,
                                      double newTime);

private:
    Constraint::Constraints& mBcs;
    std::vector<DofType> mDofs;
    std::string mSolver;
};


} // namespace NuTo
