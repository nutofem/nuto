// $Id$

#ifndef CONSTRAINTEQUATION_H_
#define CONSTRAINTEQUATION_H_

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#endif  // ENABLE_SERIALIZATION
#include <vector>

#include "nuto/mechanics/constraints/ConstraintBase.h"
#include "nuto/mechanics/constraints/ConstraintEquationTerm.h"
#include "nuto/mechanics/nodes/NodeBase.h"

namespace NuTo
{
// forward declarations
template<class T> class FullMatrix;
template<class T> class SparseMatrixCSRGeneral;

//! @brief ... constraint equations
//! @author Stefan Eckardt, ISM
//! @date 16.12.2009
class ConstraintEquation : public NuTo::ConstraintBase
{
public:
    //! @brief constructor
    //! @param rNode ... node pointer
    //! @param rDofType ... which type of dof (e.g. displacement, rotation, temperature) is constrained (first term only)
    //! @param rDofComponent ... which dof is constrained (e.g. 0 - displacement in x-direction) (first term only)
    //! @param rCoefficient ... weighting of this term in the constraint equation (first term only)
    //! @param rRhsValue ... right-hand-side value of the constraint equation
    ConstraintEquation(const NodeBase* rNode, NodeBase::eAttributes rDofType, int rDofComponent, double rCoefficient, double rRhsValue);

    //! @brief add term to constraint
    //! @param rNode ... node pointer
    //! @param rDofType ... which type of dof (e.g. displacement, rotation, temperature) is constrained
    //! @param rDofComponent ... which dof is constrained (e.g. 0 - displacement in x-direction)
    //! @param rCoefficient ... weighting of this term in the constraint equation
    void AddTerm(const NodeBase* rNode, NodeBase::eAttributes rDofType, int rDofComponent, double rCoefficient);

    //! @brief adds the constraint equation term to the matrix
    //! @param rConstraintEquation ... row in constraint matrix (is incremented during the function call)
    //! @param rConstraintMatrix ... constraint matrix
    //! @param rRHS ... right hand side vector
    void AddToConstraintMatrix(int& rConstraintEquation, NuTo::SparseMatrixCSRGeneral<double>& rConstraintMatrix, NuTo::FullMatrix<double>& rRHS) const;

    //! @brief returns the number of constraint equations
    //! @return number of constraints
    int GetNumConstraintEquations()const;

    //!@brief sets/modifies the right hand side of the constraint equations
    //!@param rRHS new right hand side
    void SetRHS(double rRHS);


#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION
private:
    std::vector<ConstraintEquationTerm> mTerms;  //!< terms of the constraint equation
    double mRhsValue;       //!< right-hand-side value of the constraint equation
};

}


#endif // CONSTRAINTEQUATION_H_
