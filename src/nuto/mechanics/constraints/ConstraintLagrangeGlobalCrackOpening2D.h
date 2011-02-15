// $Id: ConstraintLagrangeGlobalCrackOpening2D.h -1   $

#ifndef ConstraintLagrangeGlobalCrackOpening2D_H
#define ConstraintLagrangeGlobalCrackOpening2D_H

#include "nuto/mechanics/constraints/ConstraintLagrange.h"

namespace NuTo
{
template <class T>
class FullMatrix;
class StructureIp;
//! @author Joerg F. Unger, NU
//! @date June 2010
//! @brief ... class for constraints using augemented Lagrange for global crack opening to be nonnegative
class ConstraintLagrangeGlobalCrackOpening2D : public ConstraintLagrange
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
    //! @brief constructor
    //! @param rDirection ... direction of the applied constraint
    //! @param rValue ... direction of the applied constraint
    ConstraintLagrangeGlobalCrackOpening2D(const StructureIp* rStructure, double rPenaltyStiffness);

    //! @brief returns the number of constraint equations
    //! @return number of constraints
    int GetNumLagrangeMultipliers()const;

    //! @brief returns the Lagrange Multiplier
    //! first col Lagrange, second column slack variables
    void GetLagrangeMultiplier(FullMatrix<double>& rLagrangeMultiplier)const;

    //! @brief cast to Lagrange constraint - the corresponding dofs are eliminated in the global system
    NuTo::ConstraintLagrange* AsConstraintLagrange();

    //! @brief cast to Lagrange constraint - the corresponding dofs are eliminated in the global system
    const NuTo::ConstraintLagrange* AsConstraintLagrange()const;

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
    double CalculateTotalPotential()const;

    //! @brief ... print information about the object
    //! @param rVerboseLevel ... verbosity of the information
    void Info(unsigned short rVerboseLevel) const;


#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION

protected:
    //! @brief just for serialization
    ConstraintLagrangeGlobalCrackOpening2D(){};
    //! @brief structure storing the global crack
    const StructureIp *mStructure;
    //! @brief Lagrange multipliers related to the normal crack opening
    double mLagrangeValue;
    //! @brief Lagrange multipliers dofs
    int mLagrangeDOF;
};
}//namespace NuTo

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::ConstraintLagrangeGlobalCrackOpening2D)
#endif // ENABLE_SERIALIZATION

#endif //ConstraintLagrangeGlobalCrackOpening2D_H

