#include "mechanics/tools/QuasistaticProblem.h"

using namespace NuTo;

QuasistaticProblem::QuasistaticProblem(EquationSystem& s, Constraint::Constraints constraints, int numIndependentDofs,
                                       DofType dof)
    : mEquations(s)
    , mConstraints(constraints)
    , mDof(dof)
    , mCmat(mConstraints.BuildConstraintMatrix(mDof, numIndependentDofs))
{
}

void QuasistaticProblem::SetGlobalTime(double globalTime)
{
    mGlobalTime = globalTime;
}

std::pair<Eigen::SparseMatrix<double>, Eigen::VectorXd>
QuasistaticProblem::TrialSystem(const Eigen::VectorXd& x, double globalTime, double timeStep)
{
    Eigen::VectorXd deltaBrhs =
            mConstraints.GetRhs(mDof, globalTime) - mConstraints.GetRhs(mDof, globalTime + timeStep);

    auto hessian0 = mEquations.Hessian0(ToGlobalDofVector(x), {mDof});

    Eigen::VectorXd residualConstrained =
            (hessian0.JK(mDof, mDof) - mCmat.transpose() * hessian0.KK(mDof, mDof)) * deltaBrhs;

    Eigen::SparseMatrix<double> hessianMod = hessian0.JJ(mDof, mDof) - mCmat.transpose() * hessian0.KJ(mDof, mDof) -
                                             hessian0.JK(mDof, mDof) * mCmat +
                                             mCmat.transpose() * hessian0.KK(mDof, mDof) * mCmat;

    return std::make_pair(hessianMod, residualConstrained);
}

Eigen::VectorXd QuasistaticProblem::Residual(const Eigen::VectorXd& x)
{
    auto gradient = mEquations.Gradient(ToGlobalDofVector(x), {mDof});
    return gradient.J[mDof] - mCmat.transpose() * gradient.K[mDof];
}

Eigen::SparseMatrix<double> QuasistaticProblem::Derivative(const Eigen::VectorXd& x)
{
    auto hessian0 = mEquations.Hessian0(ToGlobalDofVector(x), {mDof});
    return hessian0.JJ(mDof, mDof) - mCmat.transpose() * hessian0.KJ(mDof, mDof) - hessian0.JK(mDof, mDof) * mCmat +
           mCmat.transpose() * hessian0.KK(mDof, mDof) * mCmat;
}

void QuasistaticProblem::UpdateHistory(const Eigen::VectorXd& x)
{
    mEquations.UpdateHistory(ToGlobalDofVector(x), {mDof});
}

double QuasistaticProblem::Norm(const Eigen::VectorXd& residual) const
{
    return residual.norm();
}

void QuasistaticProblem::Info(int i, const Eigen::VectorXd& x, const Eigen::VectorXd& r) const
{
    std::cout << "Iteration " << i << ": |R| = " << Norm(r) << " |x| = " << x.norm() << '\n';
}

GlobalDofVector QuasistaticProblem::ToGlobalDofVector(const Eigen::VectorXd& x)
{
    GlobalDofVector v;
    v.J[mDof] = x;
    v.K[mDof] = -mCmat * x + mConstraints.GetRhs(mDof, mGlobalTime);
    return v;
}
