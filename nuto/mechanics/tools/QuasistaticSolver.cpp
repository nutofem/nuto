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

// --------------------- ImplicitCallBack ---------------------//
ImplicitCallBack::ImplicitCallBack(TimeDependentProblem& s, ReducedSolutionSpace& reducedSolutionSpaceOperator,
                                   double tolerance)
    : mTolerance(tolerance)
    , mProblem(s)
    , mReducedSolutionSpaceOperator(reducedSolutionSpaceOperator)
{
}

Eigen::VectorXd ImplicitCallBack::Residual(const Eigen::VectorXd& u)
{
    DofVector<double> uDof;
    mReducedSolutionSpaceOperator.ToDofVector(u, uDof);

    DofVector<double> gradient =
            mProblem.Gradient(uDof, mReducedSolutionSpaceOperator.GetDofTypes(), mGlobalTime + mTimeStep, mTimeStep);
    auto f_full = ToEigen(gradient, mReducedSolutionSpaceOperator.GetDofTypes());
    return mReducedSolutionSpaceOperator.GradientToReducedBasis(f_full);
}

Eigen::SparseMatrix<double> ImplicitCallBack::Derivative(const Eigen::VectorXd& u)
{
    DofVector<double> uDof;
    mReducedSolutionSpaceOperator.ToDofVector(u, uDof);
    DofMatrixSparse<double> D_dof =
            mProblem.Hessian0(uDof, mReducedSolutionSpaceOperator.GetDofTypes(), mGlobalTime + mTimeStep, mTimeStep);
    auto D_full = ToEigen(D_dof, mReducedSolutionSpaceOperator.GetDofTypes());
    return mReducedSolutionSpaceOperator.HessianToReducedBasis(D_full);
}

void ImplicitCallBack::UpdateHistory(const Eigen::VectorXd& u)
{
    DofVector<double> uDof;
    mReducedSolutionSpaceOperator.ToDofVector(u, uDof);
    mProblem.UpdateHistory(uDof, mReducedSolutionSpaceOperator.GetDofTypes(), mGlobalTime + mTimeStep, mTimeStep);
}

double ImplicitCallBack::Norm(const Eigen::VectorXd& residual) const
{
    // all these conversions should be removed at some point, workaround to remove JK
    return residual.norm();
}

void ImplicitCallBack::Info(int i, const Eigen::VectorXd& x, const Eigen::VectorXd& r) const
{
    Log::Info << "Iteration " << i << ": |R| = " << r.norm()
              << " |x| = " << mReducedSolutionSpaceOperator.GradientToReducedBasis(x).norm() << '\n';
}

void ImplicitCallBack::SetGlobalTime(double globalTime)
{
    mTimeStep = globalTime - mGlobalTime;
    mGlobalTime = globalTime;
}

void ImplicitCallBack::SetReducedSolutionSpaceOperator(ReducedSolutionSpace& reducedSolutionSpaceOperator)
{
    mReducedSolutionSpaceOperator = reducedSolutionSpaceOperator;
}

Eigen::VectorXd ImplicitCallBack::Residual(const Eigen::VectorXd& u, double t, double dt)
{
    DofVector<double> uDof;
    mReducedSolutionSpaceOperator.ToDofVector(u, uDof);
    const std::vector<DofType>& dofTypes = mReducedSolutionSpaceOperator.GetDofTypes();
    auto gradient = mProblem.Gradient(uDof, dofTypes, t, dt);
    return ToEigen(gradient, dofTypes);
}

Eigen::SparseMatrix<double> ImplicitCallBack::Derivative(const Eigen::VectorXd& u, double t, double dt)
{
    DofVector<double> uDof;
    mReducedSolutionSpaceOperator.ToDofVector(u, uDof);
    const std::vector<DofType>& dofTypes = mReducedSolutionSpaceOperator.GetDofTypes();
    return ToEigen(mProblem.Hessian0(uDof, dofTypes, t, dt), dofTypes);
}

Eigen::VectorXd ImplicitCallBack::TrialStateRHS(const Eigen::VectorXd& u, const Eigen::SparseMatrix<double>& K_full,
                                                Eigen::VectorXd& deltaBrhsEigen, double t, double dt)
{
    DofVector<double> uDof;
    mReducedSolutionSpaceOperator.ToDofVector(u, uDof);
    Eigen::VectorXd f_full = Residual(u, t, dt);
    deltaBrhsEigen = mReducedSolutionSpaceOperator.DeltaFullRhs(mGlobalTime, t);
    return mReducedSolutionSpaceOperator.GetConstraintMatrix().transpose() * (f_full + K_full * deltaBrhsEigen);
}

// void ImplicitCallBack::FillDofVector(DofVector<double>& destination, const Eigen::VectorXd& source) const
//{
//    mReducedSolutionSpaceOperator.FillDofVector(destination, source);
//}

// --------------------- Quasistatic problem --------------------- //
void QuasistaticSolver::WriteTimeDofResidual(Eigen::VectorXd& u, std::ostream& out, DofType dofType,
                                             std::vector<int> dofNumbers, ImplicitCallBack& callBack)
{
    /* Each call to ToDofVector<double> assumes that we solve for the step t + dt. But we are in postprocess right
     * now.
     * This we are already updated to t + dt. Keeping mGlobalTime as it is, will cause the ToDofVector<double> to
     * apply
     * constraints for t + dt + dt.
     * keep in mind that the gradient is not the residual, since constraint dofs do not have a non-vanishing gradient */

    DofVector<double> uDof;
    callBack.mReducedSolutionSpaceOperator.ToDofVector(u, uDof);

    auto residual =
            callBack.mProblem.Gradient(uDof, {dofType}, callBack.mGlobalTime - callBack.mTimeStep, callBack.mTimeStep);

    double dofMean = boost::accumulate(uDof(dofType, dofNumbers), 0.) / dofNumbers.size();
    double residualSum = boost::accumulate(residual(dofType, dofNumbers), 0.);

    out << callBack.mGlobalTime << '\t' << dofMean << '\t' << residualSum << '\n';
    out << std::flush; // We really want to flush the output in case the program is interupted by CTRL-C or
    // something.
}

Eigen::VectorXd QuasistaticSolver::TrialState(Eigen::VectorXd& start, double newGlobalTime, ImplicitCallBack& callBack,
                                              std::string solverType)
{
    // update time step
    double newTimeStep = newGlobalTime - callBack.mGlobalTime;

    Eigen::VectorXd deltaBrhsEigen;

    Eigen::SparseMatrix<double> K_full = callBack.Derivative(start, callBack.mGlobalTime, callBack.mTimeStep);
    Eigen::SparseMatrix<double> Kmod = callBack.mReducedSolutionSpaceOperator.HessianToReducedBasis(K_full);
    Eigen::VectorXd fmod = callBack.TrialStateRHS(start, K_full, deltaBrhsEigen, newGlobalTime, newTimeStep);

    Eigen::VectorXd u = EigenSparseSolve(Kmod, fmod, solverType);

    u = callBack.mReducedSolutionSpaceOperator.DeltaFull(u, deltaBrhsEigen); // full

    return start - u;
}

int QuasistaticSolver::DoStep(Eigen::VectorXd& start, ImplicitCallBack& callBack, double newGlobalTime,
                              std::string solverType, double tolerance)
{
    callBack.mTolerance = tolerance; // todo: remove

    // compute trial solution (includes update of the constraint dofs, no line search)
    // trialU is with constrains
    Eigen::VectorXd trialU = TrialState(start, newGlobalTime, callBack, solverType);

    int numIterations = 0;

    // update time step
    callBack.mTimeStep = newGlobalTime - callBack.mGlobalTime;

    Eigen::VectorXd tmpX;
    try
    {
        tmpX = NewtonRaphson::Solve(
                callBack, trialU,
                NuTo::CallBackSolver(solverType, callBack.mReducedSolutionSpaceOperator.GetConstraintMatrix()), 6,
                NewtonRaphson::LineSearch(), &numIterations);
    }
    catch (std::exception& e)
    {
        throw NewtonRaphson::NoConvergence(e.what());
    }

    if (callBack.Norm(tmpX) > 1.e10)
        throw NewtonRaphson::NoConvergence("", "floating point exception");

    callBack.UpdateHistory(tmpX);
    callBack.mGlobalTime = newGlobalTime;
    start = tmpX;

    return numIterations;
}
