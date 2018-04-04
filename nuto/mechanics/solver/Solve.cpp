#include "Solve.h"
#include "nuto/math/EigenSparseSolve.h"
#include "nuto/mechanics/dofs/DofVectorConvertEigen.h"
#include "nuto/mechanics/dofs/DofMatrixSparseConvertEigen.h"

using namespace NuTo;

DofVector<double> NuTo::Solve(const DofMatrixSparse<double>& K, const DofVector<double>& f,
                              Constraint::Constraints& bcs, std::vector<DofType> dofs, std::string solver)
{
    auto K_full = ToEigen(K, dofs);
    auto f_full = ToEigen(f, dofs);

    DofMatrixSparse<double> C_dof;
    for (auto rdof : dofs)
        for (auto cdof : dofs)
            if (rdof.Id() == cdof.Id())
                C_dof(rdof, rdof) = bcs.BuildUnitConstraintMatrix2(rdof, f[rdof].rows());
            else
                C_dof(rdof, cdof) = Eigen::SparseMatrix<double>();

    auto C = ToEigen(C_dof, dofs);
    Eigen::SparseMatrix<double> Kmod = C.transpose() * K_full * C;
    Eigen::VectorXd fmod = C.transpose() * f_full;

    Eigen::VectorXd u = EigenSparseSolve(Kmod, fmod, solver);
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
    auto K_full = ToEigen(K, dofs);
    auto f_full = ToEigen(f, dofs);

    DofMatrixSparse<double> C_dof;
    for (auto rdof : dofs)
        for (auto cdof : dofs)
            if (rdof.Id() == cdof.Id())
                C_dof(rdof, rdof) = bcs.BuildUnitConstraintMatrix2(rdof, f[rdof].rows());
            else
                C_dof(rdof, cdof) = Eigen::SparseMatrix<double>();

    auto C = ToEigen(C_dof, dofs);

    // this is just for the correct size, can be replaced when the constraints know the dimensions
    DofVector<double> deltaBrhs(f);
    deltaBrhs.SetZero();
    for (auto dof : dofs)
    {
        deltaBrhs[dof] += bcs.GetSparseGlobalDeltaRhs(dof, f[dof].rows(), oldTime, newTime);
    }

    Eigen::VectorXd deltaBrhsEigen(ToEigen(deltaBrhs, dofs));
    Eigen::VectorXd residualConstrained = C.transpose() * K_full * deltaBrhsEigen;

    Eigen::SparseMatrix<double> Kmod = C.transpose() * K_full * C;
    Eigen::VectorXd fmod = C.transpose() * f_full;

    Eigen::VectorXd u = EigenSparseSolve(Kmod, fmod, solver);
    u = C * u + deltaBrhsEigen;

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

DofVector<double> ConstrainedSystemSolver::Solve(const DofMatrixSparse<double>& K, const DofVector<double>& f)
{
    return NuTo::Solve(K, f, mBcs, mDofs, mSolver);
}

DofVector<double> ConstrainedSystemSolver::SolveTrialState(const DofMatrixSparse<double>& K, const DofVector<double>& f,
                                                           double oldTime, double newTime)
{
    return NuTo::SolveTrialState(K, f, oldTime, newTime, mBcs, mDofs, mSolver);
}
