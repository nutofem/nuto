#include "nuto/mechanics/tools/QuasiStaticProblem.h"

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
QuasiStaticProblem::QuasiStaticProblem(TimeDependentProblem& s, ReducedSolutionSpace& reducedSolutionSpaceOperator,
                                       double tolerance)
    : mTolerance(tolerance)
    , mProblem(s)
    , mReducedSolutionSpaceOperator(reducedSolutionSpaceOperator)
{
}

Eigen::VectorXd QuasiStaticProblem::Residual(const Eigen::VectorXd& u)
{
    DofVector<double> uDof;
    mReducedSolutionSpaceOperator.ToDofVector(u, uDof);

    DofVector<double> gradient =
            mProblem.Gradient(uDof, mReducedSolutionSpaceOperator.GetDofTypes(), mGlobalTime + mTimeStep, mTimeStep);
    auto f_full = ToEigen(gradient, mReducedSolutionSpaceOperator.GetDofTypes());
    return mReducedSolutionSpaceOperator.GradientToReducedBasis(f_full);
}

Eigen::SparseMatrix<double> QuasiStaticProblem::Derivative(const Eigen::VectorXd& u)
{
    DofVector<double> uDof;
    mReducedSolutionSpaceOperator.ToDofVector(u, uDof);
    DofMatrixSparse<double> D_dof =
            mProblem.Hessian0(uDof, mReducedSolutionSpaceOperator.GetDofTypes(), mGlobalTime + mTimeStep, mTimeStep);
    auto D_full = ToEigen(D_dof, mReducedSolutionSpaceOperator.GetDofTypes());
    return mReducedSolutionSpaceOperator.HessianToReducedBasis(D_full);
}

void QuasiStaticProblem::UpdateHistory(const Eigen::VectorXd& u)
{
    DofVector<double> uDof;
    mReducedSolutionSpaceOperator.ToDofVector(u, uDof);
    mProblem.UpdateHistory(uDof, mReducedSolutionSpaceOperator.GetDofTypes(), mGlobalTime + mTimeStep, mTimeStep);
}

double QuasiStaticProblem::Norm(const Eigen::VectorXd& residual) const
{
    // all these conversions should be removed at some point, workaround to remove JK
    return residual.norm();
}

void QuasiStaticProblem::Info(int i, const Eigen::VectorXd& x, const Eigen::VectorXd& r) const
{
    Log::Info << "Iteration " << i << ": |R| = " << r.norm()
              << " |x| = " << mReducedSolutionSpaceOperator.GradientToReducedBasis(x).norm() << '\n';
}

void QuasiStaticProblem::SetGlobalTime(double globalTime)
{
    mTimeStep = globalTime - mGlobalTime;
    mGlobalTime = globalTime;
}

void QuasiStaticProblem::SetReducedSolutionSpaceOperator(ReducedSolutionSpace& reducedSolutionSpaceOperator)
{
    mReducedSolutionSpaceOperator = reducedSolutionSpaceOperator;
}

Eigen::VectorXd QuasiStaticProblem::Residual(const Eigen::VectorXd& u, double t, double dt)
{
    DofVector<double> uDof;
    mReducedSolutionSpaceOperator.ToDofVector(u, uDof);
    const std::vector<DofType>& dofTypes = mReducedSolutionSpaceOperator.GetDofTypes();
    auto gradient = mProblem.Gradient(uDof, dofTypes, t, dt);
    return ToEigen(gradient, dofTypes);
}

Eigen::SparseMatrix<double> QuasiStaticProblem::Derivative(const Eigen::VectorXd& u, double t, double dt)
{
    DofVector<double> uDof;
    mReducedSolutionSpaceOperator.ToDofVector(u, uDof);
    const std::vector<DofType>& dofTypes = mReducedSolutionSpaceOperator.GetDofTypes();
    return ToEigen(mProblem.Hessian0(uDof, dofTypes, t, dt), dofTypes);
}

Eigen::VectorXd QuasiStaticProblem::TrialState(Eigen::VectorXd& start, double newGlobalTime, std::string solverType)
{
    // update time step
    double newTimeStep = newGlobalTime - mGlobalTime;

    Eigen::VectorXd deltaBrhsEigen;

    Eigen::SparseMatrix<double> K_full = Derivative(start, mGlobalTime, mTimeStep);
    Eigen::VectorXd f_full = Residual(start, newGlobalTime, newTimeStep);

    deltaBrhsEigen = mReducedSolutionSpaceOperator.DeltaFullRhs(mGlobalTime, newGlobalTime);

    Eigen::SparseMatrix<double> Kmod = mReducedSolutionSpaceOperator.HessianToReducedBasis(K_full);
    Eigen::VectorXd fmod =
            mReducedSolutionSpaceOperator.GetConstraintMatrix().transpose() * (f_full + K_full * deltaBrhsEigen);

    Eigen::VectorXd u = EigenSparseSolve(Kmod, fmod, solverType);

    u = mReducedSolutionSpaceOperator.DeltaFull(u, deltaBrhsEigen); // full

    return start - u;
}
