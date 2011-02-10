// $Id$

#ifndef CONSTRAINTEQUATION_H_
#define CONSTRAINTEQUATION_H_

#include <vector>
#include "nuto/mechanics/constraints/ConstraintBase.h"
#include "nuto/mechanics/constraints/ConstraintLinear.h"
#include "nuto/mechanics/constraints/ConstraintEquationTerm.h"
#include "nuto/mechanics/nodes/NodeEnum.h"

namespace NuTo
{
// forward declarations
template<class T> class FullMatrix;
template<class T> class SparseMatrixCSRGeneral;
class NodeBase;
class ConstraintEquationTerm;
//! @brief ... constraint equations
//! @author Stefan Eckardt, ISM
//! @date 16.12.2009
class ConstraintLinearEquation :  public NuTo::ConstraintBase, public NuTo::ConstraintLinear
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION
public:
    //! @brief constructor
    //! @param rNode ... node pointer
    //! @param rDofType ... which type of dof (e.g. displacement, rotation, temperature) is constrained (first term only)
    //! @param rDofComponent ... which dof is constrained (e.g. 0 - displacement in x-direction) (first term only)
    //! @param rCoefficient ... weighting of this term in the constraint equation (first term only)
    //! @param rRhsValue ... right-hand-side value of the constraint equation
    ConstraintLinearEquation(const NodeBase* rNode, Node::eAttributes rDofType, int rDofComponent, double rCoefficient, double rRhsValue);

    //! @brief add term to constraint
    //! @param rNode ... node pointer
    //! @param rDofType ... which type of dof (e.g. displacement, rotation, temperature) is constrained
    //! @param rDofComponent ... which dof is constrained (e.g. 0 - displacement in x-direction)
    //! @param rCoefficient ... weighting of this term in the constraint equation
    void AddTerm(const NodeBase* rNode, Node::eAttributes rDofType, int rDofComponent, double rCoefficient);

    //! @brief adds the constraint equation term to the matrix
    //! @param rConstraintEquation ... row in constraint matrix (is incremented during the function call)
    //! @param rConstraintMatrix ... constraint matrix
    //! @param rRHS ... right hand side vector
    void AddToConstraintMatrix(int& rConstraintEquation, NuTo::SparseMatrixCSRGeneral<double>& rConstraintMatrix, NuTo::FullMatrix<double>& rRHS) const;

    //! @brief returns the number of constraint equations
    //! @return number of constraints
    int GetNumLinearConstraints()const;

    //! @brief cast to linear constraint - the corresponding dofs are eliminated in the global system
    NuTo::ConstraintLinear* AsConstraintLinear();

    //! @brief cast to linear constraint - the corresponding dofs are eliminated in the global system
    const NuTo::ConstraintLinear* AsConstraintLinear()const;

    //!@brief sets/modifies the right hand side of the constraint equations
    //!@param rRHS new right hand side
    void SetRHS(double rRHS);

    //! @brief ... print information about the object
    //! @param rVerboseLevel ... verbosity of the information
    void Info(unsigned short rVerboseLevel) const
    {
        throw MechanicsException("[NuTo::ConstraintLinearEquation::Info] to be implemented.");
    }

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION
protected:
    //! @brief ... just for serialize
    ConstraintLinearEquation(){};

    std::vector<ConstraintEquationTerm> mTerms;  //!< terms of the constraint equation
    double mRhsValue;       //!< right-hand-side value of the constraint equation
};

}

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::ConstraintLinearEquation)
#endif // ENABLE_SERIALIZATION

#endif // CONSTRAINTEQUATION_H_
