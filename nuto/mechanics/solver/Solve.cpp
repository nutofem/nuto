#include "Solve.h"
#include "nuto/math/EigenSparseSolve.h"

namespace NuTo
{

DofVector<double> Solve(DofMatrixSparse<double> K, DofVector<double> f, Constraint::Constraints bcs, DofType dof,
                      int numIndependentDofs, double time, std::string solver)
{
    // TODO, why do we need this method at all
    throw;

//    auto cmatUnit = bcs.BuildUnitConstraintMatrix(dof, numIndependentDofs);
//    Eigen::SparseMatrix<double> Kmod = cmatUnit.transpose() * K(dof, dof) * cmatUnit
//    Eigen::VectorXd fmod = cmatUnit.transpose() * f[dof];
//    Eigen::VectorXd fmod_constrained = (K.JK(dof, dof) - cmat.transpose() * K.KK(dof, dof)) * (-bcs.GetRhs(dof, time));

//    DofVector u;
//    u[dof] = EigenSparseSolve(Kmod, fmod + fmod_constrained, solver);
//    u.K[dof] = -cmat * u.J[dof] + bcs.GetRhs(dof, time);

//    return u;
}

} // namespace NuTo
