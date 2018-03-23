#include "nuto/mechanics/tools/QuasistaticSolver.h"

#include <ostream>
#include <boost/range/numeric.hpp>

#include "nuto/base/Timer.h"
#include "nuto/math/EigenSparseSolve.h"
#include "nuto/math/NewtonRaphson.h"
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
        mX = mProblem.RenumberDofs(constraints, mDofs, GlobalDofVector()).J;
    else
        mX = mProblem.RenumberDofs(constraints, mDofs, ToGlobalDofVector(ToEigen(mX, mDofs))).J;

    for (auto dofI : mDofs)
        for (auto dofJ : mDofs)
            if (dofI.Id() == dofJ.Id())
                mCmat(dofI, dofI) = constraints.BuildConstraintMatrix(dofI, mX[dofI].rows());
            else
                mCmat(dofI, dofJ).setZero();
}

void QuasistaticSolver::SetGlobalTime(double globalTime)
{
    mTimeStep = globalTime - mGlobalTime;
    mGlobalTime = globalTime;
}

std::pair<Eigen::SparseMatrix<double>, Eigen::VectorXd>
QuasistaticSolver::TrialSystem(const Eigen::VectorXd& x, double globalTime, double timeStep)
{
    DofVector<double> xDof = mX; // for correct size
    FromEigen(x, mDofs, &xDof);

    GlobalDofVector v;
    for (auto dof : mDofs)
    {
        v.J[dof] = xDof[dof];
        v.K[dof] = -mCmat(dof, dof) * xDof[dof] + mConstraints.GetRhs(dof, mGlobalTime);
    }

    auto hessian0 = mProblem.Hessian0(v, mDofs, globalTime, timeStep);
    auto hJJ = ToEigen(hessian0.JJ, mDofs);
    auto hJK = ToEigen(hessian0.JK, mDofs);
    auto hKJ = ToEigen(hessian0.KJ, mDofs);
    auto hKK = ToEigen(hessian0.KK, mDofs);
    auto cMat = ToEigen(mCmat, mDofs);

    NuTo::Timer myTimerStandard("standard procedure with jk components");
    Eigen::SparseMatrix<double> hessianMod = hJJ - cMat.transpose() * hKJ - hJK * cMat + cMat.transpose() * hKK * cMat;
    myTimerStandard.GetTimeDifference();
    {
        //************* check timing for computation of hessian by using H_mod = C^t H C *************

        //        //assemble hessian0Full without JK components, this if JK is completely removed, this conversion is not required, since it is
        //        //directly returned from mProblem.Hessian0
        Eigen::SparseMatrix<double> hessian0WithoutJK(hJJ.rows() + hKJ.rows(), hJJ.cols() + hJK.cols());
        hessian0WithoutJK.reserve(hJJ.nonZeros() + hJK.nonZeros() + hKJ.nonZeros() + hKK.nonZeros());
        hessian0WithoutJK.setZero();

        // Fill hessian0WithoutJK with triples from the other matrices
        std::vector<Eigen::Triplet<double> > tripletList;
        for (int k = 0; k < hJJ.outerSize(); ++k)
        {
            for (Eigen::SparseMatrix<double>::InnerIterator it(hJJ, k); it; ++it)
            {
                tripletList.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));
            }
        }
        for (int k = 0; k < hJK.outerSize(); ++k)
        {
            for (Eigen::SparseMatrix<double>::InnerIterator it(hJK, k); it; ++it)
            {
                tripletList.push_back(Eigen::Triplet<double>(it.row(), it.col()+hJJ.cols(), it.value()));
            }
        }
        for (int k = 0; k < hKJ.outerSize(); ++k)
        {
            for (Eigen::SparseMatrix<double>::InnerIterator it(hKJ, k); it; ++it)
            {
                tripletList.push_back(Eigen::Triplet<double>(it.row()+hJJ.rows(), it.col(), it.value()));
            }
        }
        for (int k = 0; k < hKK.outerSize(); ++k)
        {
            for (Eigen::SparseMatrix<double>::InnerIterator it(hKK, k); it; ++it)
            {
                tripletList.push_back(Eigen::Triplet<double>(it.row()+hJJ.rows(), it.col()+hJJ.cols(), it.value()));
            }
        }

        hessian0WithoutJK.setFromTriplets(tripletList.begin(), tripletList.end());

//        Eigen::MatrixXd hessian0WithoutJKFull(hessian0WithoutJK);
//        std::cout << "hessian0WithoutJKFull\n" << hessian0WithoutJKFull << std::endl;
//        Eigen::MatrixXd hJJFull(hJJ);
//        std::cout << "hJJFull\n" << hJJFull << std::endl;
//        Eigen::MatrixXd hJKFull(hJK);
//        std::cout << "hJKFull\n" << hJKFull << std::endl;
//        Eigen::MatrixXd hKJFull(hKJ);
//        std::cout << "hKJFull\n" << hKJFull << std::endl;
//        Eigen::MatrixXd hKKFull(hKK);
//        std::cout << "hKKFull\n" << hKKFull << std::endl;

        // compute the Cmat with the leading unit diagonal elements
        //add unit entrys, should also be done when building Cmat
        Eigen::SparseMatrix<double> cMatPlusUnit(cMat.rows()+cMat.cols(),cMat.cols());
        cMatPlusUnit.reserve(cMat.cols() + cMat.nonZeros());
        cMatPlusUnit.setZero();
        std::vector<Eigen::Triplet<double> > tripletListC;
        for (auto i=0; i<cMat.cols(); i++)
        {
            tripletListC.push_back(Eigen::Triplet<double>(i, i, 1.));
        }

        //this cMatPlusUnit = [ I]
        //                    [-C]
        for(Eigen::Index c; c<cMat.outerSize(); ++c)
        {
            for(Eigen::SparseMatrix<double>::InnerIterator itCmat(cMat, c); itCmat; ++itCmat)
                 tripletListC.push_back(Eigen::Triplet<double>(itCmat.row()+cMat.cols(), itCmat.col(), -itCmat.value()));
        }
//        cMatPlusUnit.setFromTriplets(tripletListC.begin(), tripletListC.end());

//        Eigen::MatrixXd cMatPlusUnitFull(cMatPlusUnit);
//        std::cout << "cMatPlusUnitFull\n" << cMatPlusUnitFull << std::endl;
        NuTo::Timer myTimerNew("C^T H C without jk components");
        Eigen::SparseMatrix<double> hessianModWithoutJK = cMatPlusUnit.transpose() * hessian0WithoutJK * cMatPlusUnit;
        myTimerNew.GetTimeDifference();
//        Eigen::MatrixXd hessianModWithoutJKFull(hessianModWithoutJK);
//        std::cout << "hessianModWithoutJKFull\n" << hessianModWithoutJKFull << std::endl;

//        Eigen::MatrixXd hessianModFull(hessianMod);
//        std::cout << "hessianMod\n" << hessianModFull << std::endl;

    }



    DofVector<double> deltaBrhs;
    for (auto dof : mDofs)
        deltaBrhs[dof] = mConstraints.GetRhs(dof, globalTime + timeStep) - mConstraints.GetRhs(dof, globalTime);

    Eigen::VectorXd residualConstrained = -(hJK - cMat.transpose() * hKK) * ToEigen(deltaBrhs, mDofs);

    return std::make_pair(hessianMod, residualConstrained);
}

Eigen::VectorXd QuasistaticSolver::Residual(const Eigen::VectorXd& x)
{
    auto gradient = mProblem.Gradient(ToGlobalDofVector(x), mDofs, mGlobalTime + mTimeStep, mTimeStep);
    return ToEigen(gradient.J, mDofs) - ToEigen(mCmat, mDofs).transpose() * ToEigen(gradient.K, mDofs);
}


Eigen::SparseMatrix<double> QuasistaticSolver::Derivative(const Eigen::VectorXd& x)
{
    auto hessian0 = mProblem.Hessian0(ToGlobalDofVector(x), mDofs, mGlobalTime + mTimeStep, mTimeStep);
    auto hJJ = ToEigen(hessian0.JJ, mDofs);
    auto hJK = ToEigen(hessian0.JK, mDofs);
    auto hKJ = ToEigen(hessian0.KJ, mDofs);
    auto hKK = ToEigen(hessian0.KK, mDofs);
    auto cMat = ToEigen(mCmat, mDofs);
    return hJJ - cMat.transpose() * hKJ - hJK * cMat + cMat.transpose() * hKK * cMat;
}

void QuasistaticSolver::UpdateHistory(const Eigen::VectorXd& x)
{
    mProblem.UpdateHistory(ToGlobalDofVector(x), mDofs, mGlobalTime + mTimeStep, mTimeStep);
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

GlobalDofVector QuasistaticSolver::ToGlobalDofVector(const Eigen::VectorXd& x) const
{
    DofVector<double> xDof = mX; // for correct size
    FromEigen(x, mDofs, &xDof);

    GlobalDofVector v;
    for (auto dof : mDofs)
    {
        v.J[dof] = xDof[dof];
        v.K[dof] = -mCmat(dof, dof) * xDof[dof] + mConstraints.GetRhs(dof, mGlobalTime + mTimeStep);
    }
    return v;
}

void QuasistaticSolver::WriteTimeDofResidual(std::ostream& out, DofType dofType, std::vector<int> dofNumbers)
{
    /* Disclaimer: The whole class is messy because of the time step handling via member variables. This itself is due
     * to the fact that NewtonRaphson does not care about time or time steps, but the Residual function, it tries to
     * minimize, does. So we slip that around the interface. */
    mGlobalTime -= mTimeStep;
    auto x = ToGlobalDofVector(ToEigen(mX, mDofs));
    mGlobalTime += mTimeStep;
    /* So each call to ToGlobalDofVector assumes that we solve for the step t + dt. But we are in postprocess right now.
     * This we are already updated to t + dt. Keeping mGlobalTime as it is, will cause the ToGlobalDofVector to apply
     * constraints for t + dt + dt. */

    auto residual = mProblem.Gradient(x, {dofType}, mGlobalTime, mTimeStep);

    double dofMean = boost::accumulate(x(dofType, dofNumbers), 0.) / dofNumbers.size();
    double residualSum = boost::accumulate(residual(dofType, dofNumbers), 0.);

    out << mGlobalTime << '\t' << dofMean << '\t' << residualSum << '\n';
    out << std::flush; // We really want to flush the output in case the program is interupted by CTRL-C or something.
}


int QuasistaticSolver::DoStep(double newGlobalTime, std::string solverType)
{
    EigenSparseSolver solver(solverType);

    mTimeStep = newGlobalTime - mGlobalTime;
    auto trialSystem = TrialSystem(ToEigen(mX, mDofs), mGlobalTime, mTimeStep);

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
