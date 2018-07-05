#include "Solve.h"
#include "nuto/math/EigenSparseSolve.h"
#include "nuto/mechanics/dofs/DofVectorConvertEigen.h"
#include "nuto/mechanics/dofs/DofMatrixSparseConvertEigen.h"

using namespace NuTo;

DofVector<double> NuTo::Solve(const DofMatrixSparse<double>& K, const DofVector<double>& f,
                              Constraint::Constraints& bcs, std::vector<DofType> dofs, std::string solver)
{
    auto f_full = ToEigen(f, dofs);

    DofMatrixSparse<double> C_dof;
    for (auto rdof : dofs)
        C_dof(rdof, rdof) = bcs.BuildUnitConstraintMatrix(rdof, f[rdof].rows());

    for (auto rdof : dofs)
        for (auto cdof : dofs)
            if (rdof.Id() != cdof.Id())
                C_dof(rdof, cdof) = Eigen::SparseMatrix<double>(C_dof(rdof, rdof).rows(), C_dof(cdof, cdof).cols());

    auto C = ToEigen(C_dof, dofs);
    Eigen::VectorXd fmod = C.transpose() * f_full;

    Eigen::SparseMatrix<double> Kmod;
    if (Hack::Recalculate)
    {
        auto K_full = ToEigen(K, dofs);
        Kmod = C.transpose() * K_full * C;
    }
    EigenSparseSolver s(solver);
    Eigen::VectorXd u = s.Solve(Kmod, fmod);
    u = C * u;

    // TODO: for correct size
    DofVector<double> result = f;
    FromEigen(u, f.DofTypes(), &result);
    return result;
}

DofVector<double> NuTo::SolveTrialState(const DofMatrixSparse<double>& K, const DofVector<double>& f, double oldTime,
                                        double newTime, Constraint::Constraints& bcs, std::vector<DofType> dofs,
                                        std::string solver)
{
    auto f_full = ToEigen(f, dofs);

    DofMatrixSparse<double> C_dof;
    for (auto rdof : dofs)
        C_dof(rdof, rdof) = bcs.BuildUnitConstraintMatrix(rdof, f[rdof].rows());

    for (auto rdof : dofs)
        for (auto cdof : dofs)
            if (rdof.Id() != cdof.Id())
                C_dof(rdof, cdof) = Eigen::SparseMatrix<double>(C_dof(rdof, rdof).rows(), C_dof(cdof, cdof).cols());

    auto C = ToEigen(C_dof, dofs);

    // this is just for the correct size, can be replaced when the constraints know the dimensions
    DofVector<double> deltaBrhs(f);
    deltaBrhs.SetZero();
    for (auto dof : dofs)
    {
        deltaBrhs[dof] += (bcs.GetSparseGlobalRhs(dof, f[dof].rows(), newTime) -
                           bcs.GetSparseGlobalRhs(dof, f[dof].rows(), oldTime));
    }

    Eigen::SparseMatrix<double> Kmod;
    Eigen::VectorXd fmod;
    Eigen::VectorXd deltaBrhsEigen(ToEigen(deltaBrhs, dofs));
    if (Hack::Recalculate)
    {
        auto K_full = ToEigen(K, dofs);
        Kmod = C.transpose() * K_full * C;
        fmod = C.transpose() * (f_full + K_full * deltaBrhsEigen);
    }
    else
    {
        fmod = C.transpose() * f_full;
    }

    // this last operation should in theory be done with a sparse deltaBrhsVector
    EigenSparseSolver s(solver);
    Eigen::VectorXd u = s.Solve(Kmod, fmod);
    // this is the negative increment
    // residual = gradient
    // hessian = dresidual / ddof
    // taylorseries expansion 0 = residual + hessian * deltadof
    u = C * u - deltaBrhsEigen;

    // TODO: for correct size
    DofVector<double> result = f;
    FromEigen(u, f.DofTypes(), &result);
    return result;
}

ConstrainedSystemSolver::ConstrainedSystemSolver(Constraint::Constraints& bcs, std::vector<DofType> dofs,
                                                 std::string solver)
    : mBcs(bcs)
    , mDofs(dofs)
    , mSolver(solver)
{
}

DofVector<double> ConstrainedSystemSolver::Solve(const DofMatrixSparse<double>& K, const DofVector<double>& f) const
{
    return NuTo::Solve(K, f, mBcs, mDofs, mSolver);
}

DofVector<double> ConstrainedSystemSolver::SolveTrialState(const DofMatrixSparse<double>& K, const DofVector<double>& f,
                                                           double oldTime, double newTime) const
{
    return NuTo::SolveTrialState(K, f, oldTime, newTime, mBcs, mDofs, mSolver);
}
