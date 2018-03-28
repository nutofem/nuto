#include "nuto/mechanics/tools/QuasistaticSolver.h"

#include <ostream>
#include <boost/range/numeric.hpp>

#include "nuto/base/Timer.h"
#include "nuto/math/EigenSparseSolve.h"
#include "nuto/math/NewtonRaphson.h"
#include "nuto/mechanics/dofs/DofMatrix.h"
#include "nuto/mechanics/dofs/DofVector.h"
#include "nuto/mechanics/dofs/DofVectorConvertEigen.h"
#include "nuto/mechanics/dofs/DofMatrixSparseConvertEigen.h"

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
        mX = mProblem.RenumberDofs(constraints, mDofs, DofVector<double>());
    else
        mX = mProblem.RenumberDofs(constraints, mDofs, ToDofVector(ToEigen(mX, mDofs)));

    for (auto dofI : mDofs)
        for (auto dofJ : mDofs)
            if (dofI.Id() == dofJ.Id())
                mCmatUnit(dofI, dofI) = constraints.BuildUnitConstraintMatrix2(dofI, mX[dofI].rows());
            else
                mCmatUnit(dofI, dofJ).setZero();
}

void QuasistaticSolver::SetGlobalTime(double globalTime)
{
    mTimeStep = globalTime - mGlobalTime;
    mGlobalTime = globalTime;
}

std::pair<Eigen::SparseMatrix<double>, Eigen::VectorXd> QuasistaticSolver::TrialSystem(double globalTime,
                                                                                       double timeStep)
{
    auto hessian0 = mProblem.Hessian0(mX, mDofs, globalTime, timeStep);
    auto hessian0Eigen = ToEigen(hessian0, mDofs);
    auto cMatUnit = ToEigen(mCmatUnit, mDofs);

    Eigen::SparseMatrix<double> hessianMod = cMatUnit.transpose() * hessian0Eigen * cMatUnit;

    DofVector<double> deltaBrhs;
    for (auto dof : mDofs)
    {
        deltaBrhs[dof] = mConstraints.GetSparseGlobalRhs(dof, mX[dof].rows(), globalTime + timeStep) -
                         mConstraints.GetSparseGlobalRhs(dof, mX[dof].rows(), globalTime);
    }

    Eigen::VectorXd residualConstrained = cMatUnit.transpose() * hessian0Eigen * ToEigen(deltaBrhs, mDofs);

    return std::make_pair(hessianMod, residualConstrained);
}

Eigen::VectorXd QuasistaticSolver::Residual(const Eigen::VectorXd& x)
{
    auto gradient = mProblem.Gradient(ToDofVector(x), mDofs, mGlobalTime + mTimeStep, mTimeStep);
    return ToEigen(mCmatUnit, mDofs).transpose() * ToEigen(gradient, mDofs);
}


Eigen::SparseMatrix<double> QuasistaticSolver::Derivative(const Eigen::VectorXd& x)
{
    auto hessian0 = mProblem.Hessian0(ToDofVector(x), mDofs, mGlobalTime + mTimeStep, mTimeStep);
    auto hessian0Eigen = ToEigen(hessian0, mDofs);
    auto cMatUnit = ToEigen(mCmatUnit, mDofs);
    return cMatUnit.transpose() * hessian0Eigen * cMatUnit;
}

void QuasistaticSolver::UpdateHistory(const Eigen::VectorXd& x)
{
    mProblem.UpdateHistory(ToDofVector(x), mDofs, mGlobalTime + mTimeStep, mTimeStep);
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

DofVector<double> QuasistaticSolver::ToDofVector(const Eigen::VectorXd& x) const
{
    // this result includes now all dependent and independent dofs for all active dof types
    auto cMatUnit = ToEigen(mCmatUnit, mDofs);
    Eigen::VectorXd fullX = cMatUnit * x;

    for (auto dof : mDofs)
    {
        // add the rhs of the constraint equations
        // [d_j,d_k]^T = CmatUnit * d_j (above) + b
        // since b is stored for each dof_type separately, do it in the loop
        // fullVecor.rows() is the size of the current vector including dependent and independent
        //    components of only the active dofs
        fullX += mConstraints.GetSparseGlobalRhs(dof, fullX.rows(), mGlobalTime + mTimeStep);
    }

    // for correct size, copy the inactive dofs
    DofVector<double> xDof = mX;

    // update the active dofs
    FromEigen(fullX, mDofs, &xDof);

    return xDof;
}

void QuasistaticSolver::WriteTimeDofResidual(std::ostream& out, DofType dofType, std::vector<int> dofNumbers)
{
    /* Each call to ToDofVector<double> assumes that we solve for the step t + dt. But we are in postprocess right now.
     * This we are already updated to t + dt. Keeping mGlobalTime as it is, will cause the ToDofVector<double> to apply
     * constraints for t + dt + dt. */
    auto residual = mProblem.Gradient(mX, {dofType}, mGlobalTime - mTimeStep, mTimeStep);

    double dofMean = boost::accumulate(mX(dofType, dofNumbers), 0.) / dofNumbers.size();
    double residualSum = boost::accumulate(residual(dofType, dofNumbers), 0.);

    out << mGlobalTime << '\t' << dofMean << '\t' << residualSum << '\n';
    out << std::flush; // We really want to flush the output in case the program is interupted by CTRL-C or something.
}


int QuasistaticSolver::DoStep(double newGlobalTime, std::string solverType)
{
    EigenSparseSolver solver(solverType);

    mTimeStep = newGlobalTime - mGlobalTime;
    auto trialSystem = TrialSystem(mGlobalTime, mTimeStep);

    // here is still a problem, since ToEigen(mX, mDofs) returns all dofs, whereas the solver just returns the
    // independent part
    throw;
    Eigen::VectorXd trialX = ToEigen(mX, mDofs) + solver.Solve(trialSystem.first, trialSystem.second);

    int numIterations = 0;

    Eigen::VectorXd tmpX;
    try
    {
        tmpX = NewtonRaphson::Solve(*this, trialX, solver, 6, NewtonRaphson::LineSearch(), &numIterations);
    }
    catch (std::exception& e)
    {
        throw NewtonRaphson::NoConvergence(e.what());
    }

    if (tmpX.norm() > 1.e10)
        throw NewtonRaphson::NoConvergence("", "floating point exception");

    UpdateHistory(tmpX);
    mGlobalTime = newGlobalTime;
    FromEigen(tmpX, mDofs, &mX);

    return numIterations;
}
