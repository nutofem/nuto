#include "mechanics/tools/QuasistaticSolver.h"
#include "math/EigenSparseSolve.h"
#include "math/NewtonRaphson.h"

using namespace NuTo;

QuasistaticSolver::QuasistaticSolver(TimeDependentProblem& s, DofType dof)
    : mProblem(s)
    , mDof(dof)
{
}

void QuasistaticSolver::SetConstraints(Constraint::Constraints constraints, int numIndependentDofs)
{
    mConstraints = constraints;
    mCmat = constraints.BuildConstraintMatrix(mDof, numIndependentDofs);
    mX.setZero(numIndependentDofs);
}

void QuasistaticSolver::SetGlobalTime(double globalTime)
{
    mTimeStep = globalTime - mGlobalTime;
    mGlobalTime = globalTime;
}

std::pair<Eigen::SparseMatrix<double>, Eigen::VectorXd>
QuasistaticSolver::TrialSystem(const Eigen::VectorXd& x, double globalTime, double timeStep)
{
    GlobalDofVector v;
    v.J[mDof] = x;
    v.K[mDof] = -mCmat * x + mConstraints.GetRhs(mDof, mGlobalTime);

    auto hessian0 = mProblem.Hessian0(v, {mDof}, globalTime, timeStep);
    Eigen::SparseMatrix<double> hessianMod = hessian0.JJ(mDof, mDof) - mCmat.transpose() * hessian0.KJ(mDof, mDof) -
                                             hessian0.JK(mDof, mDof) * mCmat +
                                             mCmat.transpose() * hessian0.KK(mDof, mDof) * mCmat;

    Eigen::VectorXd deltaBrhs =
            mConstraints.GetRhs(mDof, globalTime + timeStep) - mConstraints.GetRhs(mDof, globalTime);

    Eigen::VectorXd residualConstrained =
            -(hessian0.JK(mDof, mDof) - mCmat.transpose() * hessian0.KK(mDof, mDof)) * deltaBrhs;

    return std::make_pair(hessianMod, residualConstrained);
}

Eigen::VectorXd QuasistaticSolver::Residual(const Eigen::VectorXd& x)
{
    auto gradient = mProblem.Gradient(ToGlobalDofVector(x), {mDof}, mGlobalTime + mTimeStep, mTimeStep);
    return gradient.J[mDof] - mCmat.transpose() * gradient.K[mDof];
}

Eigen::SparseMatrix<double> QuasistaticSolver::Derivative(const Eigen::VectorXd& x)
{
    auto hessian0 = mProblem.Hessian0(ToGlobalDofVector(x), {mDof}, mGlobalTime + mTimeStep, mTimeStep);
    return hessian0.JJ(mDof, mDof) - mCmat.transpose() * hessian0.KJ(mDof, mDof) - hessian0.JK(mDof, mDof) * mCmat +
           mCmat.transpose() * hessian0.KK(mDof, mDof) * mCmat;
}

void QuasistaticSolver::UpdateHistory(const Eigen::VectorXd& x)
{
    mProblem.UpdateHistory(ToGlobalDofVector(x), {mDof}, mGlobalTime + mTimeStep, mTimeStep);
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
    GlobalDofVector v;
    v.J[mDof] = x;
    v.K[mDof] = -mCmat * x + mConstraints.GetRhs(mDof, mGlobalTime + mTimeStep);
    return v;
}

boost::optional<int> QuasistaticSolver::DoStep(double newGlobalTime, std::string solverType)
{
    EigenSparseSolver solver(solverType);

    mTimeStep = newGlobalTime - mGlobalTime;
    auto trialSystem = TrialSystem(mX, mGlobalTime, mTimeStep);

    Eigen::VectorXd trialX = mX + solver.Solve(trialSystem.first, trialSystem.second);

    int numIterations = 0;

    try
    {
        Eigen::VectorXd tmpX =
                NewtonRaphson::Solve(*this, trialX, solver, 12, NewtonRaphson::LineSearch(), &numIterations);
        if (tmpX.norm() > 1.e10)
            throw NewtonRaphson::NoConvergence("", "floating point exception");

        UpdateHistory(tmpX);
        mGlobalTime = newGlobalTime;
        mX = tmpX;
    }
    catch (NewtonRaphson::NoConvergence& e)
    {
        return boost::none;
    }

    return numIterations;
}
