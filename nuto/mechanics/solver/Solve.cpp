#include "Solve.h"
#include "nuto/math/EigenSparseSolve.h"
#include <Eigen/Eigenvalues>

#include <iostream>

namespace NuTo
{

GlobalDofVector Solve(GlobalDofMatrixSparse K, GlobalDofVector f, Constraint::Constraints bcs, DofType dof,
                      int numIndependentDofs, double time, std::string solver)
{
    auto cmat = bcs.BuildConstraintMatrix(dof, numIndependentDofs);
    Eigen::SparseMatrix<double> Kmod = K.JJ(dof, dof) - cmat.transpose() * K.KJ(dof, dof) - K.JK(dof, dof) * cmat +
                                       cmat.transpose() * K.KK(dof, dof) * cmat;
    Eigen::VectorXd fmod = f.J[dof] - cmat.transpose() * f.K[dof];
    Eigen::VectorXd fmod_constrained = (K.JK(dof, dof) - cmat.transpose() * K.KK(dof, dof)) * (-bcs.GetRhs(dof, time));

    GlobalDofVector u;

    Eigen::MatrixXd matrix(Kmod);
    Eigen::VectorXcd eivals = matrix.eigenvalues();
    std::cout << eivals << std::endl;

    u.J[dof] = EigenSparseSolve(Kmod, fmod + fmod_constrained, solver);
    u.K[dof] = -cmat * u.J[dof] + bcs.GetRhs(dof, time);

    return u;
}

} // namespace NuTo
