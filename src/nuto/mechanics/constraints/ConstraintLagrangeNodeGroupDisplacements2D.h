// $Id: ConstraintLagrangeNodeDisplacements2D.h -1   $

#pragma once

#include "nuto/mechanics/constraints/ConstraintLagrange.h"
#include "nuto/mechanics/constraints/ConstraintNodeGroup.h"

namespace NuTo
{
//! @author Joerg F. Unger, NU
//! @date June 2010
//! @brief ... class for all constraints applied to a node group solved using augmented Lagrange multipliers in 2D
//! additional information in: Bertsekas "constrained optimization and Lagrange multiplier methods"
class ConstraintLagrangeNodeGroupDisplacements2D : public ConstraintNodeGroup, public ConstraintLagrange
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
    //! @brief constructor
    //! @param rDirection ... direction of the applied constraint
    //! @param rValue ... direction of the applied constraint
    ConstraintLagrangeNodeGroupDisplacements2D(const Group<NodeBase>* rGroup, const NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rDirection, NuTo::Constraint::eEquationSign rEquationSign, double rValue);

    //! @brief returns the number of constraint equations
    //! @return number of constraints
    int GetNumLagrangeMultipliers()const;

    //! @brief returns the Lagrange Multiplier
    //! first col Lagrange, second column slack variables
    void GetLagrangeMultiplier(FullVector<double,Eigen::Dynamic>& rLagrangeMultiplier)const;

    //! @brief returns the Lagrange Multiplier dofs
    //! first col Lagrangedofs
    void GetDofsLagrangeMultiplier(FullVector<int,Eigen::Dynamic>& rLagrangeMultiplier)const;

    //! @brief cast to Lagrange constraint - the corresponding dofs are eliminated in the global system
    NuTo::ConstraintLagrange* AsConstraintLagrange();

    //! @brief cast to Lagrange constraint - the corresponding dofs are eliminated in the global system
    const NuTo::ConstraintLagrange* AsConstraintLagrange()const;

    //!@brief sets/modifies the right hand side of the constraint equations
    //!@param rRHS new right hand side
    void SetRHS(double rRHS);

    //! @brief sets the global dofs
    //! @param rDOF current maximum DOF, this variable is increased within the routine
    void SetGlobalDofs(int& rDOF);

    //! @brief write dof values to constraints (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    void SetGlobalDofValues(const FullVector<double,Eigen::Dynamic>& rActiveDofValues, const FullVector<double,Eigen::Dynamic>& rDependentDofValues) ;

    //! @brief extract dof values from the constraints (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    void GetGlobalDofValues(FullVector<double,Eigen::Dynamic>& rActiveDofValues, FullVector<double,Eigen::Dynamic>& rDependentDofValues) const;

    //! @brief renumber the global dofs according to predefined ordering
    //! @param rMappingInitialToNewOrdering ... mapping from initial ordering to the new ordering
    void RenumberGlobalDofs(const std::vector<int>& rMappingInitialToNewOrdering);

    //! @brief calculates the coefficient matrix for the 0-th derivative in the differential equation
    //! for a mechanical problem, this corresponds to the stiffness matrix (Hessian) of the total potential
    //! @param rResult ... coefficient matrix
    //! @param rGlobalDofs ... row and column numbers in global system
    void CalculateCoefficientMatrix_0(SparseMatrixCSRVector2Symmetric<double>& rResult,
            std::vector<int>& rGlobalDofs)const;

    //! @brief calculates the gradient of the internal potential
    //! for a mechanical problem, this corresponds to the internal force vector
    void CalculateGradientInternalPotential(NuTo::FullVector<double,Eigen::Dynamic>& rResult,
            std::vector<int>& rGlobalDofs)const;

    //! @brief calculates the internal potential
    double CalculateTotalPotential()const
    {
        throw MechanicsException("[NuTo::ConstraintLagrangeNodeGroupDisplacements2D::CalculateTotalPotential] to be implemented.");
    }

    //! @brief ... print information about the object
    //! @param rVerboseLevel ... verbosity of the information
    void Info(unsigned short rVerboseLevel) const
    {
        throw MechanicsException("[NuTo::ConstraintLagrangeNodeGroupDisplacements2D::Info] to be implemented.");
    }

    //! @brief determines the dof type affected by the constraint
    //! @return dof type
    Node::eDof GetDofType() const override
    {
        return Node::eDof::DISPLACEMENTS;
    }

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION

protected:
    //! @brief just for serialization
    ConstraintLagrangeNodeGroupDisplacements2D(){};
    //! @brief prescribed displacement of the node group (rhs)
    double mRHS;
    //! @brief direction of the applied constraint (normalized)
    double mDirection[2];

    //! @brief Lagrange multipliers related to all the nodes in the group
    std::vector<double> mLagrangeValue;
    //! @brief Lagrange multipliers dofs related to all the nodes in the group
    std::vector<double> mLagrangeDOF;
};
}//namespace NuTo

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::ConstraintLagrangeNodeGroupDisplacements2D)
#endif // ENABLE_SERIALIZATION