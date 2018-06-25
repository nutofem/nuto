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
    Eigen::VectorXd u_full = mReducedSolutionSpaceOperator.ToFull(u);
    return mReducedSolutionSpaceOperator.GradientToReducedBasis(
            ToEigen(this->FullResidual(u_full), mReducedSolutionSpaceOperator.GetDofTypes()));
}

DofVector<double> QuasiStaticProblem::FullResidual(const Eigen::VectorXd& u)
{
    mReducedSolutionSpaceOperator.ToDofVector(u, mSolution);
    return mProblem.Gradient(mSolution, mReducedSolutionSpaceOperator.GetDofTypes(), mGlobalTime + mTimeStep,
                             mTimeStep);
}

Eigen::SparseMatrix<double> QuasiStaticProblem::Derivative(const Eigen::VectorXd& u)
{
    Eigen::VectorXd u_full = mReducedSolutionSpaceOperator.ToFull(u);
    DofVector<double> uDof;
    mReducedSolutionSpaceOperator.ToDofVector(u_full, uDof);
    DofMatrixSparse<double> D_dof =
            mProblem.Hessian0(uDof, mReducedSolutionSpaceOperator.GetDofTypes(), mGlobalTime + mTimeStep, mTimeStep);
    auto D_full = ToEigen(D_dof, mReducedSolutionSpaceOperator.GetDofTypes());
    return mReducedSolutionSpaceOperator.HessianToReducedBasis(D_full);
}

void QuasiStaticProblem::UpdateHistory(const Eigen::VectorXd& u)
{
    Eigen::VectorXd u_full = mReducedSolutionSpaceOperator.ToFull(u);
    DofVector<double> uDof;
    mReducedSolutionSpaceOperator.ToDofVector(u_full, uDof);
    mProblem.UpdateHistory(uDof, mReducedSolutionSpaceOperator.GetDofTypes(), mGlobalTime + mTimeStep, mTimeStep);
}

double QuasiStaticProblem::Norm(const Eigen::VectorXd& residual) const
{
    // all these conversions should be removed at some point, workaround to remove JK
    return residual.norm();
}

void QuasiStaticProblem::Info(int i, const Eigen::VectorXd& x, const Eigen::VectorXd& r) const
{
    Log::Info << "Iteration " << i << ": |R| = " << r.norm() << " |x| = " << x.norm() << '\n';
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

void QuasiStaticProblem::FillDofVector(const Eigen::VectorXd& source, DofVector<double>& destination) const
{
    Eigen::VectorXd full = mReducedSolutionSpaceOperator.ToFull(source);
    mReducedSolutionSpaceOperator.ToDofVector(full, destination);
}

Eigen::VectorXd QuasiStaticProblem::Residual(const Eigen::VectorXd& u, double t, double dt)
{
    mReducedSolutionSpaceOperator.ToDofVector(u, mSolution);
    const std::vector<DofType>& dofTypes = mReducedSolutionSpaceOperator.GetDofTypes();
    auto gradient = mProblem.Gradient(mSolution, dofTypes, t, dt);
    return ToEigen(gradient, dofTypes);
}

Eigen::SparseMatrix<double> QuasiStaticProblem::Derivative(const Eigen::VectorXd& u, double t, double dt)
{
    mReducedSolutionSpaceOperator.ToDofVector(u, mSolution);
    const std::vector<DofType>& dofTypes = mReducedSolutionSpaceOperator.GetDofTypes();
    return ToEigen(mProblem.Hessian0(mSolution, dofTypes, t, dt), dofTypes);
}


Eigen::VectorXd QuasiStaticProblem::TrialState(DofVector<double>& start, double newGlobalTime, std::string solverType)
{
    // store old values as well as all dofs not solved for in mSolutionVector
    mSolution = start;

    // compute "old" hessian in previously equilibrated time step
    Eigen::SparseMatrix<double> K_full =
            ToEigen(mProblem.Hessian0(start, mReducedSolutionSpaceOperator.GetDofTypes(), mGlobalTime, mTimeStep),
                    mReducedSolutionSpaceOperator.GetDofTypes());

    // update time step
    mTimeStep = newGlobalTime - mGlobalTime;
    mGlobalTime = newGlobalTime;

    // compute new residual (for the new time, but with old dof values, so in most cases this is the change of external
    // force)
    Eigen::VectorXd f_full =
            ToEigen(mProblem.Gradient(start, mReducedSolutionSpaceOperator.GetDofTypes(), mGlobalTime, mTimeStep),
                    mReducedSolutionSpaceOperator.GetDofTypes());

    // compute change of rhs of constraint equations
    Eigen::VectorXd deltaBrhsEigen;
    deltaBrhsEigen = mReducedSolutionSpaceOperator.DeltaFullRhs(mGlobalTime, newGlobalTime);

    // reduce the matrix to the independent dofs
    Eigen::SparseMatrix<double> Kmod = mReducedSolutionSpaceOperator.HessianToReducedBasis(K_full);

    // compute trial resiudal
    Eigen::VectorXd fmod =
            mReducedSolutionSpaceOperator.GetConstraintMatrix().transpose() * (f_full + K_full * deltaBrhsEigen);

    // solve reduced system
    Eigen::VectorXd u = EigenSparseSolve(Kmod, fmod, solverType);

    // compute from the delta of the independent dofs the delta of dependent and independent dofs (only for active
    // doftypes)
    u = mReducedSolutionSpaceOperator.DeltaFull(u, deltaBrhsEigen);

    DofVector<double> temp;
    mReducedSolutionSpaceOperator.ToDofVector(u, temp);

    return mReducedSolutionSpaceOperator.ExtractIndependentDofVector(start - temp);
}
