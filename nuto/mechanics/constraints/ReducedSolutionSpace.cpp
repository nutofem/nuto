#include "ReducedSolutionSpace.h"
#include "nuto/mechanics/dofs/DofVectorConvertEigen.h"
#include "nuto/mechanics/dofs/DofMatrixSparseConvertEigen.h"

using namespace NuTo;
using namespace NuTo::Constraint;

ReducedSolutionSpace::ReducedSolutionSpace(const std::vector<DofType>& dofTypes, const DofContainer<int>& numTotalDofs,
                                           const Constraints& constraints)
    : mDofTypes(dofTypes)
    , mNumTotalDofs(numTotalDofs)
    , mConstraints(constraints)
{
    //        mConstraints = constraints;
    //        if (mX[mDofs.front()].rows() == 0)
    //          mX = mProblem.RenumberDofs(constraints, mDofs, DofVector<double>());
    //        else
    //            mX = mProblem.RenumberDofs(constraints, mDofs, mX);

    DofMatrixSparse<double> mCmatUnit;

    for (auto dofI : mDofTypes)
    {
        int numDofs = mNumTotalDofs.At(dofI); // throws if no dof
        for (auto dofJ : mDofTypes)
            if (dofI.Id() == dofJ.Id())
                mCmatUnit(dofI, dofI) = constraints.BuildUnitConstraintMatrix(dofI, numDofs);
            else
                mCmatUnit(dofI, dofJ).setZero();
    }

    mCmatUnitSparse = ToEigen(mCmatUnit, mDofTypes);
}

Eigen::SparseMatrix<double> ReducedSolutionSpace::HessianToReducedBasis(const Eigen::SparseMatrix<double>& matrix) const
{
    return mCmatUnitSparse.transpose() * matrix * mCmatUnitSparse;
}

Eigen::VectorXd ReducedSolutionSpace::GradientToReducedBasis(Eigen::VectorXd& gradient) const
{
    return mCmatUnitSparse.transpose() * gradient;
}

Eigen::VectorXd ReducedSolutionSpace::ToFull(const Eigen::VectorXd& independent) const
{
    return mCmatUnitSparse * independent;
}

Eigen::VectorXd ReducedSolutionSpace::ToFullWithRhs(const Eigen::VectorXd& independent, double time) const
{
    DofType dof = mDofTypes[0]; // see constructor, so far only one entry
    return mCmatUnitSparse * independent + mConstraints.GetRhs(dof, time);
}

DofVector<double> ReducedSolutionSpace::ToDofVector(const Eigen::VectorXd& full,
                                                    const DofVector<double>& dofVector) const
{
    DofVector<double> dofVectorNew(dofVector);
    dofVectorNew[mDofTypes[0]] = full;
}

Eigen::VectorXd ReducedSolutionSpace::DeltaFull(const Eigen::VectorXd& independent,
                                                const Eigen::VectorXd& deltaBrhsEigen) const
{
    return mCmatUnitSparse * independent - deltaBrhsEigen;
}

Eigen::VectorXd ReducedSolutionSpace::DeltaFullRhs(double timeOld, double timeNew) const
{
    DofVector<double> deltaBrhs;
    for (auto dof : mDofTypes)
    {
        deltaBrhs[dof].setZero(mNumTotalDofs[dof]);
        deltaBrhs[dof] += (mConstraints.GetSparseGlobalRhs(dof, mNumTotalDofs[dof], timeNew) -
                           mConstraints.GetSparseGlobalRhs(dof, mNumTotalDofs[dof], timeOld));
    }

    return ToEigen(deltaBrhs, mDofTypes);
}

Eigen::VectorXd ReducedSolutionSpace::ToReducedBasis(const DofVector<double>& dofVector, double time) const
{
    return mCmatUnitSparse * dofVector[mDofTypes[0]];
}
