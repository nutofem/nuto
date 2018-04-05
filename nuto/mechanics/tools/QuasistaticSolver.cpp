#include "nuto/mechanics/tools/QuasistaticSolver.h"

#include <ostream>
#include <iomanip>
#include <boost/range/numeric.hpp>

#include "nuto/base/Timer.h"
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
        mX = mProblem.RenumberDofs(constraints, mDofs, mX);

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

DofVector<double> QuasistaticSolver::TrialState(double newGlobalTime, const ConstrainedSystemSolver& solver)
{
    // compute hessian for last converged time step
    auto hessian0 = mProblem.Hessian0(mX, mDofs, mGlobalTime, mTimeStep);
    Eigen::MatrixXd hessian0Eigen(ToEigen(hessian0, mDofs));

    // update time step
    double newTimeStep = newGlobalTime - mGlobalTime;

    // compute residual for new time step (in particular, these are changing external forces)
    auto gradient = mProblem.Gradient(mX, mDofs, newGlobalTime, newTimeStep);

    DofVector<double> trialU = mX - solver.SolveTrialState(hessian0, gradient, mGlobalTime, newGlobalTime);

    return trialU;
}

DofVector<double> QuasistaticSolver::Residual(const DofVector<double>& u)
{
    return mProblem.Gradient(u, mDofs, mGlobalTime + mTimeStep, mTimeStep);
}


DofMatrixSparse<double> QuasistaticSolver::Derivative(const DofVector<double>& u)
{
    return mProblem.Hessian0(u, mDofs, mGlobalTime + mTimeStep, mTimeStep);
}

void QuasistaticSolver::UpdateHistory(const DofVector<double>& x)
{
    mProblem.UpdateHistory(x, mDofs, mGlobalTime + mTimeStep, mTimeStep);
}

double QuasistaticSolver::Norm(const DofVector<double>& residual) const
{
    // all these conversions should be removed at some point, workaround to remove JK
    Eigen::VectorXd tmp = ToEigen(residual, mDofs);
    auto C = ToEigen(mCmatUnit, mDofs);
    return (C.transpose() * tmp).norm();
}

void QuasistaticSolver::Info(int i, const DofVector<double>& x, const DofVector<double>& r) const
{
    if (mQuiet)
        return;
    std::cout << std::right << std::setfill(' ');
    std::cout << "Iteration " << i << ": |R| = " << std::setw(11) << Norm(r) << " |x| = " << std::setw(11) << Norm(x)
              << '\n';
}

void QuasistaticSolver::WriteTimeDofResidual(std::ostream& out, DofType dofType, std::vector<int> dofNumbers)
{
    /* Each call to ToDofVector<double> assumes that we solve for the step t + dt. But we are in postprocess right
     * now.
     * This we are already updated to t + dt. Keeping mGlobalTime as it is, will cause the ToDofVector<double> to
     * apply
     * constraints for t + dt + dt.
     * keep in mind that the gradient is not the residual, since constraint dofs do not have a non-vanishing gradient
*/
    auto residual = mProblem.Gradient(mX, {dofType}, mGlobalTime - mTimeStep, mTimeStep);

    double dofMean = boost::accumulate(mX(dofType, dofNumbers), 0.) / dofNumbers.size();
    double residualSum = boost::accumulate(residual(dofType, dofNumbers), 0.);

    out << mGlobalTime << '\t' << dofMean << '\t' << residualSum << '\n';
    out << std::flush; // We really want to flush the output in case the program is interupted by CTRL-C or
    // something.
}

int QuasistaticSolver::DoStep(double newGlobalTime, std::string solverType)
{
    // allocate constraint system solver
    ConstrainedSystemSolver solver(mConstraints, mDofs, solverType);

    // compute trial solution (includes update of the constraint dofs, no line search)
    DofVector<double> trialU = TrialState(newGlobalTime, solver);

    int numIterations = 0;

    // update time step
    mTimeStep = newGlobalTime - mGlobalTime;

    DofVector<double> tmpX;
    try
    {
        tmpX = NewtonRaphson::Solve(*this, trialU, solver, 6, NewtonRaphson::LineSearch(), &numIterations);
    }
    catch (std::exception& e)
    {
        throw NewtonRaphson::NoConvergence(e.what());
    }

    if (Norm(tmpX) > 1.e10)
        throw NewtonRaphson::NoConvergence("", "floating point exception");

    UpdateHistory(tmpX);
    mGlobalTime = newGlobalTime;
    mX = tmpX;

    return numIterations;
}
