#pragma once

#include "nuto/base/Exception.h"
#include "nuto/mechanics/constraints/Constraints.h"

namespace NuTo
{

class LinearConstraintApplicator
{
public:
    LinearConstraintApplicator(const std::vector<DofType>& dofTypes, const DofContainer<int>& mNumTotalDofs,
                               const Constraints& constraints)
        : mDofTypes(dofTypes)
        , mDofContainer<int> mNumTotalDofs
        , mConstraints(constraints)
    {
        //this is a temporary fix, usually the CMat does not have to be build doftype-wise
        DofMatrixSparse<double> C_dof;

        for (auto dofI : mDofTypes)
            for (auto dofJ : mDofTypes)
                if (dofI.Id() == dofJ.Id())
                    C_dof(dofI, dofI) = constraints.BuildUnitConstraintMatrix(dofI, mNumTotalDofs[dofI]);
                else
                    C_dof(dofI, dofJ).setZero();

        mCmat = ToEigen(C_dof, mDofTypes);
    }

    //! @brief transforms a hessian matrix from the full system by computing C^T K C to the reduced system
    //! @param Mat matrix in the original space
    //! @return hessian in the reduced space
    Eigen::SparseMatrix<double> HessianToReducedBasis(const Eigen::SparseMatrix<double>& m) override
    {
        return mCmat.transpose() * m * mCmat;
    }

    //! @brief transforms a gradient from the full system by computing C^T v to the reduced system
    //! @param v vector in the original space
    //! @return hessian in the reduced space
    Eigen::SparseMatrix<double> GradientToReducedBasis(const Eigen::SparseMatrix<double>& Mat) override
    Eigen::VectorXd GradientToReducedBasis(Eigen::VectorXd& v) override
    {
        return mCmat.transpose()*v;
    }

    //! @brief transforms a reduced solution vector back to the full system
    //! @param uj vector of independent dofs
    //! @param t time to compute the dependent dofs based on the constraint equations
    //! @return C*uj+b(t)
    Eigen::VectorXd ToFull(const Eigen::VectorXd& uj, double t) override
    {
        return mCmat*uj+b(t);
    }

    //! @brief transforms a delta of the reduced solution vector back to the full system (assuming t=const)
    //! @param uj vector of independent dofs
    //! @return C*duj
    Eigen::VectorXd ToDeltaFull(const Eigen::VectorXd& uj) override
    {
        return mCmat*uj;
    }

private:
    // these mdofTypes and the mNumTotalDofs are only tmp, there should be a reference to the function space that is
    // constraint by the LinearConstraintapplicator
    //! @brief doftypes used in the solution procedure
    std::vector<DofType> mDofTypes;
    //! @brief number of total dofs for each dof type in the solution procedure
    DofContainer<int> mNumTotalDofs;
    //! @brief the linear transformation matrix to map the reduced space to the full space d_tot = C d_red + b(t)
    Eigen::Sparse<double> mCmat;
    //! @brief constraints
    Constraints mConstraints;
};

} // namespace NuTo
