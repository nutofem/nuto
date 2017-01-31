#include <iostream>

#include "mechanics/MechanicsException.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/groups/Group.h"
#include "mechanics/constraints/ConstraintLinearNodeGroupRotations2D.h"
#include "math/SparseMatrixCSRGeneral.h"

//! @brief constructor
NuTo::ConstraintLinearNodeGroupRotations2D::ConstraintLinearNodeGroupRotations2D(const Group<NodeBase>* rGroup, double rValue) :
        ConstraintNodeGroup(rGroup), ConstraintLinear()
{
    mRHS = rValue;
}

//! @brief returns the number of constraint equations
//! @return number of constraints
int NuTo::ConstraintLinearNodeGroupRotations2D::GetNumLinearConstraints()const
{
    return mGroup->GetNumMembers();
}

//!@brief sets/modifies the right hand side of the constraint equation
//!@param rRHS new right hand side
void NuTo::ConstraintLinearNodeGroupRotations2D::SetRHS(double rRHS)
{
	mRHS=rRHS;
}

//! @brief adds the constraint equations to the matrix
//! @param curConstraintEquation (is incremented during the function call)
//! @param rConstraintMatrix (the first row where a constraint equation is added is given by curConstraintEquation)
void NuTo::ConstraintLinearNodeGroupRotations2D::AddToConstraintMatrix(int& curConstraintEquation,
        NuTo::SparseMatrix<double>& rConstraintMatrix)const
{
    for (Group<NodeBase>::const_iterator itNode=mGroup->begin(); itNode!=mGroup->end(); itNode++)
    {
        if (itNode->second->GetNum(Node::eDof::ROTATIONS)!=1)
            throw MechanicsException("[NuTo::ConstraintLinearNodeGroupRotations2D::AddToConstraintMatrix] Node should have exactly 1 rotational dof.");

        rConstraintMatrix.AddValue(curConstraintEquation,itNode->second->GetDof(Node::eDof::ROTATIONS, 0),1);

        curConstraintEquation++;
    }
}

//!@brief writes for the current constraint equation(s) the rhs into the vector
// (in case of more than one equation per constraint, curConstraintEquation is increased based on the number of constraint equations per constraint)
//! @param curConstraintEquation (is incremented during the function call)
//! @param rConstraintMatrix (the first row where a constraint equation is added is given by curConstraintEquation)
void NuTo::ConstraintLinearNodeGroupRotations2D::GetRHS(int& curConstraintEquation,Eigen::VectorXd& rRHS)const
{
    for (Group<NodeBase>::const_iterator itNode=mGroup->begin(); itNode!=mGroup->end(); itNode++)
    {
        rRHS(curConstraintEquation,0) = mRHS;
        if (itNode->second->GetNum(Node::eDof::ROTATIONS)!=1)
            throw MechanicsException("[NuTo::ConstraintLinearNodeGroupRotations2D::GetRHS] Node should have exactly 1 rotational dof.");

        curConstraintEquation++;
    }
}

NuTo::Node::eDof NuTo::ConstraintLinearNodeGroupRotations2D::GetDofType() const
{
    return Node::eDof::ROTATIONS;
}

#ifdef ENABLE_SERIALIZATION
// serialize
template<class Archive>
void NuTo::ConstraintLinearNodeGroupRotations2D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ConstraintLinearNodeGroupRotations2D" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstraintNodeGroup)
       & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstraintLinear)
       & BOOST_SERIALIZATION_NVP(mRHS);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ConstraintLinearNodeGroupRotations2D" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::ConstraintLinearNodeGroupRotations2D)
#endif // ENABLE_SERIALIZATION