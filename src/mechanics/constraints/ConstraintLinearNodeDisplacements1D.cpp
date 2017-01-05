// $Id$
#include <iostream>

#include "mechanics/MechanicsException.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/constraints/ConstraintLinearNodeDisplacements1D.h"
#include "math/FullMatrix.h"
#include "math/SparseMatrixCSRGeneral.h"

// constructor
NuTo::ConstraintLinearNodeDisplacements1D::ConstraintLinearNodeDisplacements1D(const NodeBase* rNode, double rDirection, double rValue):
        ConstraintNode(rNode), ConstraintLinear()
{
    // set direction
    if (std::abs(rDirection) < 1e-14)
    {
        throw MechanicsException("[NuTo::ConstraintLinearNodeDisplacements1D::ConstraintLinearNodeDisplacements1D] Length of the direction vector is zero");
    }
    // set normalized direction
    this->mDirection = rDirection / std::abs(rDirection);

    // set value
    mRHS = rValue;
}

//! @brief returns the number of constraint equations
//! @return number of constraints
int NuTo::ConstraintLinearNodeDisplacements1D::GetNumLinearConstraints()const
{
    return 1;
}

//!@brief sets/modifies the right hand side of the constraint equation
//!@param rRHS new right hand side
void NuTo::ConstraintLinearNodeDisplacements1D::SetRHS(double rRHS)
{
	mRHS=rRHS;
}

//! @brief adds the constraint equations to the matrix
//! @param curConstraintEquation (is incremented during the function call)
//! @param rConstraintMatrix (the first row where a constraint equation is added is given by curConstraintEquation)
void NuTo::ConstraintLinearNodeDisplacements1D::AddToConstraintMatrix(int& curConstraintEquation,
        NuTo::SparseMatrix<double>& rConstraintMatrix)const
{
    // add constraint to constrain matrix
    if (mNode->GetNum(Node::eDof::DISPLACEMENTS)!=1)
    {
        throw MechanicsException("[NuTo::ConstraintLinearNodeDisplacements1D::ConstraintBase] Node does not have displacements or has more than one displacement component.");
    }
    rConstraintMatrix.AddValue(curConstraintEquation,mNode->GetDof(Node::eDof::DISPLACEMENTS, 0), this->mDirection);

    // increase constraint equation number
    curConstraintEquation++;
}

//!@brief writes for the current constraint equation(s) the rhs into the vector
// (in case of more than one equation per constraint, curConstraintEquation is increased based on the number of constraint equations per constraint)
//! @param curConstraintEquation (is incremented during the function call)
//! @param rConstraintMatrix (the first row where a constraint equation is added is given by curConstraintEquation)
void NuTo::ConstraintLinearNodeDisplacements1D::GetRHS(int& curConstraintEquation,Eigen::VectorXd& rRHS)const
{
    // add constraint to constrain matrix
    if (mNode->GetNum(Node::eDof::DISPLACEMENTS)!=1)
    {
        throw MechanicsException("[NuTo::ConstraintLinearNodeDisplacements1D::ConstraintBase] Node does not have displacements or has more than one displacement component.");
    }
    // set right hand side value
    rRHS(curConstraintEquation) = mRHS;

    // increase constraint equation number
    curConstraintEquation++;
}

NuTo::Node::eDof NuTo::ConstraintLinearNodeDisplacements1D::GetDofType() const
{
    return Node::eDof::DISPLACEMENTS;
}


#ifdef ENABLE_SERIALIZATION
// serialize
template<class Archive>
void NuTo::ConstraintLinearNodeDisplacements1D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ConstraintLinearNodeDisplacements1D" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstraintNode);
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstraintLinear);
    ar & BOOST_SERIALIZATION_NVP(mRHS);
    ar & BOOST_SERIALIZATION_NVP(mDirection);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ConstraintLinearNodeDisplacements1D" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::ConstraintLinearNodeDisplacements1D)

void NuTo::ConstraintLinearNodeDisplacements1D::SetNodePtrAfterSerialization(const std::map<uintptr_t, uintptr_t>& mNodeMapCast)
{
    NuTo::ConstraintNode::SetNodePtrAfterSerialization(mNodeMapCast);
}

#endif // ENABLE_SERIALIZATION
