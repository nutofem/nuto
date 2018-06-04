#pragma once
#include <eigen3/Eigen/Sparse>
#include "nuto/mechanics/dofs/DofContainer.h"
#include "nuto/mechanics/dofs/DofVector.h"
#include "nuto/mechanics/dofs/DofMatrixSparse.h"
#include "Constraints.h"

namespace NuTo
{
namespace Constraint
{
class ReducedSolutionSpace
{
public:
    ReducedSolutionSpace()
    {
    }

    ReducedSolutionSpace(const std::vector<DofType>& dofTypes, const DofContainer<int>& numTotalDofs,
                         const Constraints& constraints);

    //! @brief transforms a matrix from the full system to the reduced system
    //! @param matrix ... in the original space
    //! @return C^T * matrix * C
    Eigen::SparseMatrix<double> HessianToReducedBasis(const Eigen::SparseMatrix<double>& matrix) const;

    //! @brief transforms a gradient from the full system by computing C^T v to the reduced system
    //! @param gradient ... vector in the original space
    //! @return C^T * gradient
    Eigen::VectorXd GradientToReducedBasis(const Eigen::VectorXd& gradient) const;

    //! @brief transforms a reduced solution vector back to the full system
    //! @param independent ... vector of independent dofs
    //! @return C*independent
    Eigen::VectorXd ToFull(const Eigen::VectorXd& independent) const;

    //! @brief transforms a reduced solution vector back to the full system with the time dependent constraints
    //! @param independent ... vector of independent dofs
    //! @param time ... to compute the dependent dofs based on the constraint equations
    //! @return C*independent + RHS(time)
    Eigen::VectorXd ToFullWithRhs(const Eigen::VectorXd& independent, double time) const;

    //! @brief transforms a full solution vector back to a dofVector and replaces the corresponding entries
    //! @param full ... vector of all active dofs
    //! @param dofVector ... dof vector of all dofs
    void ToDofVector(const Eigen::VectorXd& source, DofVector<double>& destination) const;

    //! @brief transforms a delta of the reduced solution vector back to the full system (assuming t=const)
    //! @param independent ... vector of independent dofs
    //! @return C*independent
    Eigen::VectorXd DeltaFull(const Eigen::VectorXd& independent, const Eigen::VectorXd& deltaBrhsEigen) const;

    Eigen::VectorXd DeltaFullRhs(double timeOld, double timeNew) const;

    void FillDofVector(DofVector<double>& destination, const Eigen::VectorXd& source) const;

    //! @brief reduces the dof solution vector to the vector of independent dofs
    //! @param dofVector ... actual state of the system
    //! @param time ... in order to verify consistent initial conditions (constraint equations are fulfilled)
    //! @return vector ... with only the independent dofs
    Eigen::VectorXd ToReducedBasis(const DofVector<double>& dofVector) const;

    const std::vector<DofType>& GetDofTypes() const
    {
        return mDofTypes;
    }

    const Eigen::SparseMatrix<double>& GetConstraintMatrix() const
    {
        return mCmatUnitSparse;
    }

private:
    // these mdofTypes and the mNumTotalDofs are only tmp, there should be a reference to the function space that is
    // constraint by the ReducedSolutionSpace

    //! @brief doftypes used in the solution procedure
    std::vector<DofType> mDofTypes;

    //! @brief number of total dofs for each dof type in the solution procedure
    DofContainer<int> mNumTotalDofs;

    //! @brief the linear transformation matrix to map the reduced space to the full space d_tot = C d_red + b(t)
    Eigen::SparseMatrix<double> mCmatUnitSparse;

    DofMatrixSparse<double> mCmatUnit;

    //! @brief constraints
    Constraints mConstraints;
};
} // -- NuTo
} // -- Constraint
