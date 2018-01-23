#include "Solve.h"
#include "EigenSolver.h"

bool isPrefix(std::string input, std::string prefix)
{
    return input.compare(0, prefix.length(), prefix) == 0;
}


namespace NuTo
{

Eigen::VectorXd Solve(Eigen::SparseMatrix<double> A, Eigen::VectorXd b, std::string solver)
{
    if (isPrefix(solver, "Eigen"))
        return EigenSolver(A, b, solver);
    throw Exception("Que?");
}


GlobalDofVector Solve(GlobalDofMatrixSparse K, GlobalDofVector f, Constraint::Constraints bcs, DofType dof,
                      int numIndependentDofs, double time, std::string solver)
{
    auto cmat = bcs.BuildConstraintMatrix(dof, numIndependentDofs);
    Eigen::SparseMatrix<double> Kmod = K.JJ(dof, dof) - cmat.transpose() * K.KJ(dof, dof) - K.JK(dof, dof) * cmat +
                                       cmat.transpose() * K.KK(dof, dof) * cmat;
    Eigen::VectorXd fmod = f.J[dof] - cmat.transpose() * f.K[dof];
    Eigen::VectorXd fmod_constrained = (K.JK(dof, dof) - cmat.transpose() * K.KK(dof, dof)) * (-bcs.GetRhs(dof, time));

    GlobalDofVector u;
    u.J[dof] = Solve(Kmod, fmod + fmod_constrained, solver);
    u.K[dof] = -cmat * u.J[dof] + bcs.GetRhs(dof, time);

    return u;
}

} // namespace NuTo
