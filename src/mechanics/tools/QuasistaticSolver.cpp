#include "mechanics/tools/QuasistaticSolver.h"
#include "math/EigenSparseSolve.h"
#include "math/NewtonRaphson.h"
#include "mechanics/dofs/DofVectorConvertEigen.h"
#include "mechanics/dofs/DofMatrixSparseConvertEigen.h"

using namespace NuTo;

QuasistaticSolver::QuasistaticSolver(TimeDependentProblem& s, DofType dof)
    : mProblem(s)
    , mDofs({dof})
{
}

QuasistaticSolver::QuasistaticSolver(TimeDependentProblem& s, std::vector<DofType> dofs)
    : mProblem(s)
    , mDofs(dofs)
{
}

void QuasistaticSolver::SetConstraints(Constraint::Constraints constraints)
{
    mConstraints = constraints;
    if (mX[mDofs.front()].rows() == 0)
        mX = mProblem.RenumberDofs(constraints, mDofs, GlobalDofVector()).J;
    else
        mX = mProblem.RenumberDofs(constraints, mDofs, ToGlobalDofVector(ToEigen(mX, mDofs))).J;

    for (auto dofI : mDofs)
        for (auto dofJ : mDofs)
            if (dofI.Id() == dofJ.Id())
                mCmat(dofI, dofI) = constraints.BuildConstraintMatrix(dofI, mX[dofI].rows());
            else
                mCmat(dofI, dofJ).setZero();
}

void QuasistaticSolver::SetGlobalTime(double globalTime)
{
    mTimeStep = globalTime - mGlobalTime;
    mGlobalTime = globalTime;
}

std::pair<Eigen::SparseMatrix<double>, Eigen::VectorXd>
QuasistaticSolver::TrialSystem(const Eigen::VectorXd& x, double globalTime, double timeStep)
{
    DofVector<double> xDof = mX; // for correct size
    FromEigen(x, mDofs, &xDof);

    GlobalDofVector v;
    for (auto dof : mDofs)
    {
        v.J[dof] = xDof[dof];
        v.K[dof] = -mCmat(dof, dof) * xDof[dof] + mConstraints.GetRhs(dof, mGlobalTime);
    }

    auto hessian0 = mProblem.Hessian0(v, mDofs, globalTime, timeStep);
    auto hJJ = ToEigen(hessian0.JJ, mDofs);
    auto hJK = ToEigen(hessian0.JK, mDofs);
    auto hKJ = ToEigen(hessian0.KJ, mDofs);
    auto hKK = ToEigen(hessian0.KK, mDofs);
    auto cMat = ToEigen(mCmat, mDofs);
    Eigen::SparseMatrix<double> hessianMod = hJJ - cMat.transpose() * hKJ - hJK * cMat + cMat.transpose() * hKK * cMat;

    DofVector<double> deltaBrhs;
    for (auto dof : mDofs)
        deltaBrhs[dof] = mConstraints.GetRhs(dof, globalTime + timeStep) - mConstraints.GetRhs(dof, globalTime);

    Eigen::VectorXd residualConstrained = -(hJK - cMat.transpose() * hKK) * ToEigen(deltaBrhs, mDofs);

    return std::make_pair(hessianMod, residualConstrained);
}

Eigen::VectorXd QuasistaticSolver::Residual(const Eigen::VectorXd& x)
{
    auto gradient = mProblem.Gradient(ToGlobalDofVector(x), mDofs, mGlobalTime + mTimeStep, mTimeStep);
    return ToEigen(gradient.J, mDofs) - ToEigen(mCmat, mDofs).transpose() * ToEigen(gradient.K, mDofs);
}

Eigen::SparseMatrix<double> QuasistaticSolver::Derivative(const Eigen::VectorXd& x)
{
    auto hessian0 = mProblem.Hessian0(ToGlobalDofVector(x), mDofs, mGlobalTime + mTimeStep, mTimeStep);
    auto hJJ = ToEigen(hessian0.JJ, mDofs);
    auto hJK = ToEigen(hessian0.JK, mDofs);
    auto hKJ = ToEigen(hessian0.KJ, mDofs);
    auto hKK = ToEigen(hessian0.KK, mDofs);
    auto cMat = ToEigen(mCmat, mDofs);
    return hJJ - cMat.transpose() * hKJ - hJK * cMat + cMat.transpose() * hKK * cMat;
}

void QuasistaticSolver::UpdateHistory(const Eigen::VectorXd& x)
{
    mProblem.UpdateHistory(ToGlobalDofVector(x), mDofs, mGlobalTime + mTimeStep, mTimeStep);
}

double QuasistaticSolver::Norm(const Eigen::VectorXd& residual) const
{
    return residual.cwiseAbs().maxCoeff();
}

void QuasistaticSolver::Info(int i, const Eigen::VectorXd& x, const Eigen::VectorXd& r) const
{
    if (mQuiet)
        return;
    std::cout << "Iteration " << i << ": |R| = " << Norm(r) << " |x| = " << x.norm() << '\n';
}

GlobalDofVector QuasistaticSolver::ToGlobalDofVector(const Eigen::VectorXd& x)
{
    DofVector<double> xDof = mX; // for correct size
    FromEigen(x, mDofs, &xDof);

    GlobalDofVector v;
    for (auto dof : mDofs)
    {
        v.J[dof] = xDof[dof];
        v.K[dof] = -mCmat(dof, dof) * xDof[dof] + mConstraints.GetRhs(dof, mGlobalTime + mTimeStep);
    }
    return v;
}

int QuasistaticSolver::DoStep(double newGlobalTime, std::string solverType)
{
    EigenSparseSolver solver(solverType);

    mTimeStep = newGlobalTime - mGlobalTime;
    auto trialSystem = TrialSystem(ToEigen(mX, mDofs), mGlobalTime, mTimeStep);

    Eigen::VectorXd trialX = ToEigen(mX, mDofs) + solver.Solve(trialSystem.first, trialSystem.second);

    int numIterations = 0;

    Eigen::VectorXd tmpX = NewtonRaphson::Solve(*this, trialX, solver, 12, NewtonRaphson::LineSearch(), &numIterations);
    if (tmpX.norm() > 1.e10)
        throw NewtonRaphson::NoConvergence("", "floating point exception");

    UpdateHistory(tmpX);
    mGlobalTime = newGlobalTime;
    FromEigen(tmpX, mDofs, &mX);

    return numIterations;
}
