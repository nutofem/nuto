#include "nuto/mechanics/tools/NewmarkSolver.h"

#include "nuto/math/EigenSparseSolve.h"
#include "nuto/mechanics/dofs/DofMatrixSparseConvertEigen.h"
#include "nuto/mechanics/dofs/DofVectorConvertEigen.h"

namespace NuTo
{

NewmarkSolver::NewmarkSolver(TimeDependentProblem<2>& equations, NuTo::DofType dof, double beta, double gamma)
    : mProblem{equations}
    , mDofs{{dof}}
    , mBeta{beta}
    , mGamma{gamma}
{
}

NewmarkSolver::NewmarkSolver(TimeDependentProblem<2>& equations, std::vector<NuTo::DofType> dofs, double beta,
                             double gamma)
    : mProblem{equations}
    , mDofs{dofs}
    , mBeta{beta}
    , mGamma{gamma}
{
}

void NewmarkSolver::SetConstraints(NuTo::Constraint::Constraints constraints)
{
    mConstraints = constraints;
}

void NewmarkSolver::SetGlobalTime(double globalTime)
{
    mTimeStep = globalTime - mGlobalTime;
    mGlobalTime = globalTime;
}

Eigen::SparseMatrix<double> NewmarkSolver::Hessian(const Eigen::SparseMatrix<double>& hessian0,
                                                   const Eigen::SparseMatrix<double>& hessian1,
                                                   const Eigen::SparseMatrix<double>& hessian2, double delta_t) const
{
    return hessian0 + mGamma / (delta_t * mBeta) * hessian1 + 1 / (delta_t * delta_t * mBeta) * hessian2;
}

DofVector<double> NewmarkSolver::TrialState(double newGlobalTime, const NuTo::ConstrainedSystemSolver& solver)
{
    // update time step
    double delta_t = newGlobalTime - mGlobalTime;

    // compute hessians for last converged time step
    Eigen::SparseMatrix<double> hessian0 = ToEigen(mProblem.Hessian<0>(mX, mDofs, mGlobalTime, mTimeStep), mDofs);
    Eigen::SparseMatrix<double> hessian1 = ToEigen(mProblem.Hessian<1>(mDofs, mGlobalTime, mTimeStep), mDofs);
    Eigen::SparseMatrix<double> hessian2 = ToEigen(mProblem.Hessian<2>(mDofs, mGlobalTime, mTimeStep), mDofs);

    // compute effective hessian
    Eigen::SparseMatrix<double> hessian = Hessian(hessian0, hessian1, hessian2, delta_t);


    // compute residual for new time step (in particular, these are changing external forces)
    DofVector<double> gradient_dof = mProblem.Gradient(mDofs, newGlobalTime, delta_t);
    Eigen::VectorXd gradient = ToEigen(gradient_dof, mDofs);


    //    auto K_full = ToEigen(K, dofs);
    //    auto f_full = ToEigen(f, dofs);

    DofMatrixSparse<double> C_dof;
    for (auto rdof : mDofs)
        C_dof(rdof, rdof) = mConstraints.BuildUnitConstraintMatrix(rdof, gradient_dof[rdof].rows());

    for (auto rdof : mDofs)
        for (auto cdof : mDofs)
            if (rdof.Id() != cdof.Id())
                C_dof(rdof, cdof) = Eigen::SparseMatrix<double>(C_dof(rdof, rdof).rows(), C_dof(cdof, cdof).cols());

    auto C = ToEigen(C_dof, mDofs);

    // this is just for the correct size, can be replaced when the constraints know the dimensions
    DofVector<double> deltaBrhs_dof(gradient_dof);
    deltaBrhs_dof.SetZero();
    for (auto dof : mDofs)
    {
        deltaBrhs_dof[dof] += (mConstraints.GetSparseGlobalRhs(dof, gradient_dof[dof].rows(), newGlobalTime) -
                               mConstraints.GetSparseGlobalRhs(dof, gradient_dof[dof].rows(), mGlobalTime));
    }

    Eigen::SparseMatrix<double> hessianMod = C.transpose() * hessian * C;

    Eigen::VectorXd deltaBrhs(ToEigen(deltaBrhs_dof, mDofs));

    // this last operation should in theory be done with a sparse deltaBrhsVector

    Eigen::VectorXd v = ToEigen(mX[1], mDofs);
    Eigen::VectorXd a = ToEigen(mX[2], mDofs);

    Eigen::VectorXd fMod = C.transpose() * (gradient + hessian * deltaBrhs -
                                            hessian1 * (mGamma / mBeta * v - delta_t * (1 - mGamma / (2 * mBeta)) * a) -
                                            hessian2 * (1 / (delta_t * mBeta) * v - 1 / (2 * mBeta) * a));

    // !!! Fixed Solver set here. This should be corrected later. Problem: The passed solver takes DofVectors/Matrices
    // and does additional stuff instead of just solving. There is no Getter for the internally used solver string.
    Eigen::VectorXd u = EigenSparseSolve(hessianMod, fMod, "EigenSparseLU");

    // this is the negative increment
    // residual = gradient
    // hessian = dresidual / ddof
    // taylorseries expansion 0 = residual + hessian * deltadof
    u = C * u - deltaBrhs;

    // TODO: for correct size
    DofVector<double> result = gradient_dof;
    FromEigen(u, result.DofTypes(), &result);
    //    return result;


    DofVector<double> trialU = mX[0] - result;

    return trialU;
}

void NewmarkSolver::UpdateHistory(const DofVector<double>& x)
{
    mProblem.UpdateHistory(x, mDofs, mGlobalTime + mTimeStep, mTimeStep);
}

DofVector<double> NewmarkSolver::Residual(const DofVector<double>& u)
{
    return mProblem.Gradient(u, mDofs, mGlobalTime + mTimeStep, mTimeStep);
}

DofMatrixSparse<double> NewmarkSolver::Derivative(const DofVector<double>& u)
{
    Eigen::SparseMatrix<double> hessian0 = ToEigen(mProblem.Hessian<0>(mX, mDofs, mGlobalTime, mTimeStep), mDofs);
    Eigen::SparseMatrix<double> hessian1 = ToEigen(mProblem.Hessian<1>(mDofs, mGlobalTime, mTimeStep), mDofs);
    Eigen::SparseMatrix<double> hessian2 = ToEigen(mProblem.Hessian<2>(mDofs, mGlobalTime, mTimeStep), mDofs);

    DofMatrixSparse<double> hessian_dof;
    // FromEigen(Hessian(hessian0, hessian1, hessian2, mTimeStep), mDofs, &hessian_dof); //! <--- No such function -.-
    throw Exception(__PRETTY_FUNCTION__, "Implementation missing.");
    return hessian_dof;
}


} // namespace NuTo
