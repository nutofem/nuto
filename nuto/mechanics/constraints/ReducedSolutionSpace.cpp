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
    DofMatrixSparse<double> CmatUnit;
    DofMatrixSparse<double> CmatUnitInv;
    for (auto dofI : mDofTypes)
    {
        int numDofsI = mNumTotalDofs.At(dofI); // throws if no dof
        CmatUnit(dofI, dofI) = constraints.BuildUnitConstraintMatrix(dofI, numDofsI);
        CmatUnitInv(dofI, dofI) = constraints.BuildUnitConstraintMatrixInv(dofI, numDofsI);
    }

    for (auto dofI : mDofTypes)
    {
        int numDofsI = mNumTotalDofs.At(dofI); // throws if no dof
        for (auto dofJ : mDofTypes)
        {
            if (dofI.Id() != dofJ.Id())
            {
                CmatUnit(dofI, dofJ) =
                        Eigen::SparseMatrix<double>(CmatUnit(dofI, dofI).rows(), CmatUnit(dofJ, dofJ).cols());
                CmatUnitInv(dofI, dofJ) =
                        Eigen::SparseMatrix<double>(CmatUnit(dofJ, dofJ).cols(), CmatUnit(dofI, dofI).rows());
            }
        }
    }

    mCmatUnitSparse = ToEigen(CmatUnit, mDofTypes);
    mCmatUnitSparseInv = ToEigen(CmatUnitInv, mDofTypes);
}


Eigen::SparseMatrix<double> ReducedSolutionSpace::HessianToReducedBasis(const Eigen::SparseMatrix<double>& matrix) const
{
    return mCmatUnitSparse.transpose() * matrix * mCmatUnitSparse;
}

Eigen::VectorXd ReducedSolutionSpace::GradientToReducedBasis(const Eigen::VectorXd& gradient) const
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

void ReducedSolutionSpace::ToDofVector(const Eigen::VectorXd& source, DofVector<double>& destination) const
{
    for (auto dof : mDofTypes)
        destination[dof].setZero(mNumTotalDofs[dof]);

    FromEigen(source, mDofTypes, &destination);
}

Eigen::VectorXd ReducedSolutionSpace::ExtractIndependentDofVector(const DofVector<double>& source) const
{
    return mCmatUnitSparseInv * ToEigen(source, mDofTypes);
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

void ReducedSolutionSpace::FillDofVector(DofVector<double>& destination, const Eigen::VectorXd& source) const
{
    for (auto dof : mDofTypes)
        destination[dof].setZero(mNumTotalDofs[dof]);

    FromEigen(source, mDofTypes, &destination);
}


Eigen::VectorXd ReducedSolutionSpace::ToReducedBasis(const DofVector<double>& dofVector) const
{
    Eigen::VectorXd full = ToEigen(dofVector, mDofTypes);
    return mCmatUnitSparse.transpose() * full;
}
