#include "nuto/mechanics/tools/QuasistaticSolver.h"

#include <ostream>
#include <iomanip>
#include <boost/range/numeric.hpp>

#include "nuto/base/Logger.h"
#include "nuto/math/NewtonRaphson.h"
#include "nuto/mechanics/dofs/DofMatrix.h"
#include "nuto/mechanics/dofs/DofVector.h"
#include "nuto/mechanics/dofs/DofVectorConvertEigen.h"
#include "nuto/mechanics/dofs/DofMatrixSparseConvertEigen.h"

using namespace NuTo;

// --------------------- NewtonCallBack ---------------------//
NewtonCallBack::NewtonCallBack(TimeDependentProblem& s, ReducedSolutionSpace& reducedSolutionSpaceOperator,
                               double tolerance)
    : mTolerance(tolerance)
    , mProblem(s)
    , mReducedSolutionSpaceOperator(reducedSolutionSpaceOperator)
{
}

DofVector<double> NewtonCallBack::Residual(double globalTime, double timeStep)
{
    return mProblem.Gradient(mReducedSolutionSpaceOperator.GetDofTypes(), globalTime + timeStep, timeStep);
}


DofMatrixSparse<double> NewtonCallBack::Derivative(double globalTime, double timeStep)
{
    return mProblem.Hessian0(mReducedSolutionSpaceOperator.GetDofTypes(), globalTime + timeStep, timeStep);
}

void NewtonCallBack::Update(const DofVector<double>& state, double globalTime, double timeStep)
{
    mProblem.Update(state, mReducedSolutionSpaceOperator.GetDofTypes(), globalTime + timeStep, timeStep);
}

double NewtonCallBack::Norm(const DofVector<double>& residual) const
{
    // all these conversions should be removed at some point, workaround to remove JK
    Eigen::VectorXd tmp = ToEigen(residual, mReducedSolutionSpaceOperator.GetDofTypes());
    auto C = mReducedSolutionSpaceOperator.GetConstraintMatrix();
    return (C.transpose() * tmp).norm();
}

void NewtonCallBack::Info(int i, const DofVector<double>& x, const DofVector<double>& r) const
{
    Log::Info << "Iteration " << i << ": |R| = " << Norm(r) << " |x| = " << Norm(x) << '\n';
}

// --------------------- Quasistatic problem --------------------- //
void QuasistaticSolver::SetGlobalTime(double globalTime)
{
    mTimeStep = globalTime - mGlobalTime;
    mGlobalTime = globalTime;
}

void QuasistaticSolver::WriteTimeDofResidual(std::ostream& out, DofType dofType, std::vector<int> dofNumbers,
                                             TimeDependentProblem& problem)
{
    /* Each call to ToDofVector<double> assumes that we solve for the step t + dt. But we are in postprocess right
     * now.
     * This we are already updated to t + dt. Keeping mGlobalTime as it is, will cause the ToDofVector<double> to
     * apply
     * constraints for t + dt + dt.
     * keep in mind that the gradient is not the residual, since constraint dofs do not have a non-vanishing gradient
*/
    auto residual = problem.Gradient(problem.GetDofState(), {dofType}, mGlobalTime - mTimeStep, mTimeStep);

    double dofMean = boost::accumulate(problem.GetDofState()(dofType, dofNumbers), 0.) / dofNumbers.size();
    double residualSum = boost::accumulate(residual(dofType, dofNumbers), 0.);

    out << mGlobalTime << '\t' << dofMean << '\t' << residualSum << '\n';
    out << std::flush; // We really want to flush the output in case the program is interupted by CTRL-C or
    // something.
}

DofVector<double> QuasistaticSolver::TrialState(double newGlobalTime, NewtonCallBack& problem, std::string solverType)
{
    // compute hessian for last converged time step
    auto hessian0 = problem.Derivative(mGlobalTime, mTimeStep);

    // update time step
    double newTimeStep = newGlobalTime - mGlobalTime;

    // compute residual for new time step (in particular, these are changing external forces)
    auto gradient = problem.Residual(newGlobalTime, newTimeStep);

    Eigen::SparseMatrix<double> K_full = ToEigen(hessian0, problem.GetReducedSolutionSpaceOperator().GetDofTypes());
    Eigen::VectorXd f_full = ToEigen(gradient, problem.GetReducedSolutionSpaceOperator().GetDofTypes());

    Eigen::SparseMatrix<double> Kmod = problem.GetReducedSolutionSpaceOperator().HessianToReducedBasis(K_full);
    Eigen::VectorXd deltaBrhsEigen = problem.GetReducedSolutionSpaceOperator().DeltaFullRhs(mGlobalTime, newGlobalTime);
    Eigen::VectorXd fmod = problem.GetReducedSolutionSpaceOperator().GetConstraintMatrix().transpose() *
                           (f_full + K_full * deltaBrhsEigen);

    Eigen::VectorXd u = EigenSparseSolve(Kmod, fmod, solverType);

    u = problem.GetReducedSolutionSpaceOperator().DeltaFull(u, deltaBrhsEigen);

    // TODO: for correct size
    DofVector<double> result = gradient;
    FromEigen(u, gradient.DofTypes(), &result);

    DofVector<double> trialU = problem.GetProblem().GetDofState() - result;

    return trialU;
}

int QuasistaticSolver::DoStep(TimeDependentProblem& problem, ReducedSolutionSpace& reducedSolutionSpaceOperator,
                              double newGlobalTime, std::string solverType, double tolerance)
{
    NewtonCallBack newtonCallBack(problem, reducedSolutionSpaceOperator, tolerance);

    // compute trial solution (includes update of the constraint dofs, no line search)
    DofVector<double> trialU = TrialState(newGlobalTime, newtonCallBack, solverType);

    int numIterations = 0;

    // update time step
    mTimeStep = newGlobalTime - mGlobalTime;

    DofVector<double> tmpX;
    try
    {
        tmpX = NewtonRaphson::Solve(newtonCallBack, trialU, mGlobalTime + mTimeStep, mTimeStep, solverType, 6,
                                    NewtonRaphson::LineSearch(), &numIterations);
    }
    catch (std::exception& e)
    {
        throw NewtonRaphson::NoConvergence(e.what());
    }

    if (newtonCallBack.Norm(tmpX) > 1.e10)
        throw NewtonRaphson::NoConvergence("", "floating point exception");

    newtonCallBack.Update(tmpX, mGlobalTime + mTimeStep, mTimeStep);
    mGlobalTime = newGlobalTime;

    return numIterations;
}
