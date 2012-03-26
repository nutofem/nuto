// $Id$

#ifndef CONSTRAINTEQUATIONTERM_H_
#define CONSTRAINTEQUATIONTERM_H_

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif  // ENABLE_SERIALIZATION

#include "nuto/mechanics/nodes/NodeEnum.h"

namespace NuTo
{
class NodeBase;
template<class T> class SparseMatrixCSRGeneral;

//! @brief ... term in a constraint equation
//! @author Stefan Eckardt, ISM
//! @date 16.12.2009
class ConstraintEquationTerm
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION
public:
    //! @brief constructor
    //! @param rNode ... node pointer
    //! @param rDofType ... which type of dof (e.g. displacement, rotation, temperature) is constrained
    //! @param rDofComponent ... which dof is constrained (e.g. 0 - dispalacement in x-direction)
    //! @param rCoefficient ... weighting of this term in the constraint equation
    ConstraintEquationTerm(const NodeBase* rNode, Node::eAttributes rDofType, int rDofComponent, double rCoefficient);

    //! @brief returns the dof that is related to that linear term
    //! @return dof
    int GetDof() const;

    //! @brief returns the coefficient that is related to that linear term
    //! @return coefficient
    double GetCoefficient() const
    {
    	return mCoefficient;
    }

    //! @brief adds the constraint equation term to the matrix
    //! @param rRow ... row in constraint matrix
    //! @param rConstraintMatrix ... constraint matrix
    void AddToConstraintMatrix(int rRow, NuTo::SparseMatrixCSRGeneral<double>& rConstraintMatrix) const;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION
protected:
    const NodeBase* mNode;           //!< node pointer
    Node::eAttributes mDofType;  //!< which type of dof (e.g. displacement, rotation, temperature) is constrained
    int mDofComponent;               //!< which dof is constrained (e.g. 0 - dispalacement in x-direction)
    double mCoefficient;             //!< weighting of this term in the constraint equation

    //! @brief default constructor
    ConstraintEquationTerm();
};

}

#endif // CONSTRAINTEQUATIONTERM_H_
