#include "mechanics/tools/QuasistaticProblem.h"

using namespace NuTo;

QuasistaticProblem::QuasistaticProblem(TimeDependentProblem& s, DofType dof)
    : mProblem(s)
    , mDof(dof)
{
}

void QuasistaticProblem::SetConstraints(Constraint::Constraints constraints, int numIndependentDofs)
{
    mConstraints = constraints;
    mCmat = constraints.BuildConstraintMatrix(mDof, numIndependentDofs);
}

void QuasistaticProblem::SetGlobalTime(double globalTime)
{
    mGlobalTime = globalTime;
}

std::pair<Eigen::SparseMatrix<double>, Eigen::VectorXd>
QuasistaticProblem::TrialSystem(const Eigen::VectorXd& x, double globalTime, double timeStep)
{
    Eigen::VectorXd deltaBrhs =
            mConstraints.GetRhs(mDof, globalTime + timeStep) - mConstraints.GetRhs(mDof, globalTime);

    auto hessian0 = mProblem.Hessian0(ToGlobalDofVector(x), {mDof}, globalTime, timeStep);

    Eigen::VectorXd residualConstrained =
            -(hessian0.JK(mDof, mDof) - mCmat.transpose() * hessian0.KK(mDof, mDof)) * deltaBrhs;

    Eigen::SparseMatrix<double> hessianMod = hessian0.JJ(mDof, mDof) - mCmat.transpose() * hessian0.KJ(mDof, mDof) -
                                             hessian0.JK(mDof, mDof) * mCmat +
                                             mCmat.transpose() * hessian0.KK(mDof, mDof) * mCmat;

    return std::make_pair(hessianMod, residualConstrained);
}

Eigen::VectorXd QuasistaticProblem::Residual(const Eigen::VectorXd& x)
{
    auto gradient = mProblem.Gradient(ToGlobalDofVector(x), {mDof}, mGlobalTime, 0.);
    return gradient.J[mDof] - mCmat.transpose() * gradient.K[mDof];
}

Eigen::SparseMatrix<double> QuasistaticProblem::Derivative(const Eigen::VectorXd& x)
{
    auto hessian0 = mProblem.Hessian0(ToGlobalDofVector(x), {mDof}, mGlobalTime, 0.);
    return hessian0.JJ(mDof, mDof) - mCmat.transpose() * hessian0.KJ(mDof, mDof) - hessian0.JK(mDof, mDof) * mCmat +
           mCmat.transpose() * hessian0.KK(mDof, mDof) * mCmat;
}

void QuasistaticProblem::UpdateHistory(const Eigen::VectorXd& x)
{
    mProblem.UpdateHistory(ToGlobalDofVector(x), {mDof}, mGlobalTime, 0.);
}

double QuasistaticProblem::Norm(const Eigen::VectorXd& residual) const
{
    return residual.cwiseAbs().maxCoeff();
}

void QuasistaticProblem::Info(int i, const Eigen::VectorXd& x, const Eigen::VectorXd& r) const
{
    if (mQuiet)
        return;
    std::cout << "Iteration " << i << ": |R| = " << Norm(r) << " |x| = " << x.norm() << '\n';
}

GlobalDofVector QuasistaticProblem::ToGlobalDofVector(const Eigen::VectorXd& x)
{
    GlobalDofVector v;
    v.J[mDof] = x;
    v.K[mDof] = -mCmat * x + mConstraints.GetRhs(mDof, mGlobalTime);
    return v;
}
