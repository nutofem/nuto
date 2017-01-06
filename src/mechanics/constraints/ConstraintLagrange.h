// $Id: ConstraintLagrange.h 328 2010-10-01 14:39:32Z unger3 $

#pragma once

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif  // ENABLE_SERIALIZATION

#include <vector>
#include "mechanics/constraints/ConstraintNonlinear.h"

namespace NuTo
{
class ConstraintLinear;
template<class T> class SparseMatrixCSRVector2Symmetric;

namespace Constraint
{
    enum class eEquationSign;
}//namespace Constraint

//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... standard abstract class for all constraint equations solved with augmented Lagrange
class ConstraintLagrange : public ConstraintNonlinear
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
    //! @brief constructor
    ConstraintLagrange(NuTo::Constraint::eEquationSign rEquationSign);

    //! @brief constructor
    ConstraintLagrange(NuTo::Constraint::eEquationSign rEquationSign, double rPenaltyStiffness);

    //! @brief destructor
    virtual ~ConstraintLagrange();

    //! @brief returns the Lagrange Multiplier
    //! first col Lagrange
    virtual void GetLagrangeMultiplier(Eigen::VectorXd& rLagrangeMultiplier)const=0;

    //! @brief returns the Lagrange Multiplier dofs
    //! first col Lagrangedofs
    virtual void GetDofsLagrangeMultiplier(Eigen::VectorXd& rLagrangeMultiplier)const=0;

    //! @brief sets the global dofs
    //! @param rDOF current maximum DOF, this variable is increased within the routine
    virtual void SetGlobalDofs(int& rDOF)=0;

    //! @brief sets the penalty stiffness of the augmented lagrangian
    //! @param rDOF current maximum DOF, this variable is increased within the routine
    void SetPenaltyStiffness(double rPenalty);

    //! @brief write dof values to constraints (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void SetGlobalDofValues(const Eigen::VectorXd& rActiveDofValues, const Eigen::VectorXd& rDependentDofValues) = 0;

    //! @brief extract dof values from the constraints (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void GetGlobalDofValues(Eigen::VectorXd& rActiveDofValues, Eigen::VectorXd& rDependentDofValues) const = 0;

    //! @brief renumber the global dofs according to predefined ordering
    //! @param rMappingInitialToNewOrdering ... mapping from initial ordering to the new ordering
    virtual void RenumberGlobalDofs(const std::vector<int>& rMappingInitialToNewOrdering) = 0;

    //! @brief calculates the coefficient matrix for the 0-th derivative in the differential equation
    //! for a mechanical problem, this corresponds to the stiffness matrix (Hessian) of the total potential
    //! @param rResult ... coefficient matrix
    //! @param rGlobalDofs ... row and column numbers in global system
    virtual void CalculateCoefficientMatrix_0(NuTo::SparseMatrixCSRVector2Symmetric<double>& rResult,
            std::vector<int>& rGlobalDofs)const=0;

    //! @brief calculates the gradient of the internal potential
    //! for a mechanical problem, this corresponds to the internal force vector
    virtual void CalculateGradientInternalPotential(Eigen::VectorXd& rResult,
            std::vector<int>& rGlobalDofs)const=0;

    //! @brief calculates the internal potential
    virtual double CalculateTotalPotential()const=0;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION


protected:
    // just for serialization
    ConstraintLagrange(){};
    //for the Lagrange method, also inequations can be defined, '<', '=' or '>'
    Constraint::eEquationSign mEquationSign;
    //penalty stiffness for augmented Lagrange
    double mPenalty;

};
}//namespace NuTo
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::ConstraintLagrange)
#endif // ENABLE_SERIALIZATION


