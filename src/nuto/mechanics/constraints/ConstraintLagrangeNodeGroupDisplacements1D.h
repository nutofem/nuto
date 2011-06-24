// $Id: ConstraintLagrangeNodeDisplacements2D.h -1   $

#ifndef CONSTRAINTNODEDISPLACEMENTS1D_H
#define CONSTRAINTLAGRANGENODEDISPLACEMENTS1D_H

#include "nuto/mechanics/constraints/ConstraintLagrange.h"
#include "nuto/mechanics/constraints/ConstraintNodeGroup.h"

namespace NuTo
{
template <class T>
class FullMatrix;
//! @author Joerg F. Unger, NU
//! @date June 2010
//! @brief ... class for all constraints applied to a node group solved using Lagrange multipliers in 1D
//! additional information in: Bertsekas "constrained optimization and Lagrange multiplier methods"
class ConstraintLagrangeNodeGroupDisplacements1D : public ConstraintNodeGroup, public ConstraintLagrange
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
    //! @brief constructor
    //! @param rDirection ... direction of the applied constraint
    //! @param rValue ... direction of the applied constraint
    ConstraintLagrangeNodeGroupDisplacements1D(const Group<NodeBase>* rGroup, const NuTo::FullMatrix<double>& rDirection, NuTo::Constraint::eEquationSign rEquationSign, double rValue);

    //! @brief returns the number of constraint equations
    //! @return number of constraints
    int GetNumLagrangeMultipliers()const;

    //! @brief returns the Lagrange Multiplier
    //! first col Lagrange, second column slack variables
    void GetLagrangeMultiplier(FullMatrix<double>& rLagrangeMultiplier)const;

    //! @brief returns the Lagrange Multiplier dofs
    //! first col Lagrangedofs
    void GetDofsLagrangeMultiplier(FullMatrix<int>& rLagrangeMultiplier)const;

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
    void SetGlobalDofValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues) ;

    //! @brief extract dof values from the constraints (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    void GetGlobalDofValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const;

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
    void CalculateGradientInternalPotential(NuTo::FullMatrix<double>& rResult,
            std::vector<int>& rGlobalDofs)const;

    //! @brief calculates the internal potential
    double CalculateTotalPotential()const
    {
        throw MechanicsException("[NuTo::ConstraintLagrangeNodeGroupDisplacements1D::CalculateTotalPotential] to be implemented.");
    }

    //! @brief ... print information about the object
    //! @param rVerboseLevel ... verbosity of the information
    void Info(unsigned short rVerboseLevel) const
    {
        throw MechanicsException("[NuTo::ConstraintLagrangeNodeGroupDisplacements1D::Info] to be implemented.");
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
    ConstraintLagrangeNodeGroupDisplacements1D(){};
    //! @brief prescribed displacement of the node group (rhs)
    double mRHS;

    //! @brief Lagrange multipliers related to all the nodes in the group
    std::vector<double> mLagrangeValue;
    //! @brief Lagrange multipliers dofs related to all the nodes in the group
    std::vector<double> mLagrangeDOF;
};
}//namespace NuTo

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::ConstraintLagrangeNodeGroupDisplacements1D)
#endif // ENABLE_SERIALIZATION

#endif //CONSTRAINTLAGRANGENODEDISPLACEMENTS1D_H

