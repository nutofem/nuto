#pragma once
#include <Eigen/Dense>
#include <Eigen/SparseLU>
#include "mechanics/dofs/GlobalDofVector.h"
#include "mechanics/dofs/GlobalDofMatrixSparse.h"

namespace NuTo
{

Eigen::VectorXd Solve(Eigen::SparseMatrix<double> A, Eigen::VectorXd b)
{
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.analyzePattern(A);
    solver.factorize(A);
    return solver.solve(b); 
}

// TODO: DofType, numIndependentDofs and time need to go eventually
GlobalDofVector Solve(GlobalDofMatrixSparse K, GlobalDofVector f, Constraint::Constraints bcs, DofType dof,
                      int numIndependentDofs, double time)
{
    auto cmat = bcs.BuildConstraintMatrix(dof, numIndependentDofs);
    Eigen::SparseMatrix<double> Kmod = K.JJ(dof, dof) - cmat.transpose() * K.KJ(dof, dof) - K.JK(dof, dof) * cmat +
                                       cmat.transpose() * K.KK(dof, dof) * cmat;
    Eigen::VectorXd fmod = f.J[dof] - cmat.transpose() * f.K[dof];
    Eigen::VectorXd fmod_constrained = (K.JK(dof, dof) - cmat.transpose() * K.KK(dof, dof)) * (-bcs.GetRhs(dof, time));

    GlobalDofVector u;
    u.J[dof] = Solve(Kmod, fmod + fmod_constrained);
    u.K[dof] = -cmat * u.J[dof] + bcs.GetRhs(dof, time);

    return u;
}

} // namespace NuTo
