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
template<class T> class SparseMatrixCSRGeneral;
class NodeBase;
class ConstraintEquationTerm;
//! @brief ... constraint equations
//! @author Stefan Eckardt, ISM
//! @date 16.12.2009
class ConstraintLinearEquation :  public NuTo::ConstraintLinear
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION
    friend class NewmarkIndirect;
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
    void AddToConstraintMatrix(int& rConstraintEquation, NuTo::SparseMatrixCSRGeneral<double>& rConstraintMatrix) const;

    //!@brief writes for the current constraint equation(s) the rhs into the vector
    // (in case of more than one equation per constraint, curConstraintEquation is increased based on the number of constraint equations per constraint)
    //! @param curConstraintEquation (is incremented during the function call)
    //! @param rConstraintMatrix (the first row where a constraint equation is added is given by curConstraintEquation)
    void GetRHS(int& curConstraintEquation,NuTo::FullVector<double,Eigen::Dynamic>& rRHS)const;

    //!@brief returns the rhs
    double GetRHS()const;

    //! @brief returns the number of constraint equations
    //! @return number of constraints
    int GetNumLinearConstraints()const;

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
