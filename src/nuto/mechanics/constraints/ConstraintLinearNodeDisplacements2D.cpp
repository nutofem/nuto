// $Id$

#include <iostream>

#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/constraints/ConstraintLinearNodeDisplacements2D.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"

NuTo::ConstraintLinearNodeDisplacements2D::ConstraintLinearNodeDisplacements2D(const NodeBase* rNode, const NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rDirection, double rValue) :
        ConstraintNode(rNode), ConstraintLinear()
{
    if (rDirection.GetNumColumns()!=1 || rDirection.GetNumRows()!=2)
        throw MechanicsException("[NuTo::ConstraintLinearNodeDisplacements2D::ConstraintLinearNodeDisplacements2D] Dimension of the direction matrix must be equal to the dimension of the structure.");

    memcpy(mDirection,rDirection.data(),2*sizeof(double));
    //normalize the direction
    double norm = sqrt(mDirection[0]*mDirection[0]+mDirection[1]*mDirection[1]);
    if (norm < 1e-14)
    {
        throw MechanicsException("[NuTo::ConstraintLinearNodeDisplacements2D::ConstraintLinearNodeDisplacements2D] direction vector has zero length.");
    }
    double invNorm = 1./norm;
    mDirection[0]*=invNorm;
    mDirection[1]*=invNorm;
    mRHS = rValue;
}

//! @brief returns the number of constraint equations
//! @return number of constraints
int NuTo::ConstraintLinearNodeDisplacements2D::GetNumLinearConstraints()const
{
    return 1;
}

//!@brief sets/modifies the right hand side of the constraint equation
//!@param rRHS new right hand side
void NuTo::ConstraintLinearNodeDisplacements2D::SetRHS(double rRHS)
{
	mRHS=rRHS;
}

//! @brief adds the constraint equations to the matrix
//! @param curConstraintEquation (is incremented during the function call)
//! @param rConstraintMatrix (the first row where a constraint equation is added is given by curConstraintEquation)
void NuTo::ConstraintLinearNodeDisplacements2D::AddToConstraintMatrix(int& curConstraintEquation,
        NuTo::SparseMatrix<double>& rConstraintMatrix)const
{
    if (mNode->GetNum(Node::DISPLACEMENTS)!=2)
        throw MechanicsException("[NuTo::ConstraintLinearNodeDisplacements2D::ConstraintBase] Node does not have displacements or has more than one displacement component.");
    if (std::abs(mDirection[0])>1e-18)
        rConstraintMatrix.AddValue(curConstraintEquation,mNode->GetDof(Node::DISPLACEMENTS, 0),mDirection[0]);
    if (std::abs(mDirection[1])>1e-18)
        rConstraintMatrix.AddValue(curConstraintEquation,mNode->GetDof(Node::DISPLACEMENTS, 1),mDirection[1]);

    curConstraintEquation++;
}

//!@brief writes for the current constraint equation(s) the rhs into the vector
// (in case of more than one equation per constraint, curConstraintEquation is increased based on the number of constraint equations per constraint)
//! @param curConstraintEquation (is incremented during the function call)
//! @param rConstraintMatrix (the first row where a constraint equation is added is given by curConstraintEquation)
void NuTo::ConstraintLinearNodeDisplacements2D::GetRHS(int& curConstraintEquation,NuTo::FullVector<double,Eigen::Dynamic>& rRHS)const
{
    if (mNode->GetNum(Node::DISPLACEMENTS)!=2)
        throw MechanicsException("[NuTo::ConstraintLinearNodeDisplacements2D::ConstraintBase] Node does not have displacements or has more than one displacement component.");

    rRHS(curConstraintEquation) = mRHS;

    curConstraintEquation++;
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template<class Archive>
void NuTo::ConstraintLinearNodeDisplacements2D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ConstraintLinearNodeDisplacements2D" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstraintNode);
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstraintLinear);
    ar & BOOST_SERIALIZATION_NVP(mRHS);
    ar & BOOST_SERIALIZATION_NVP(mDirection);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ConstraintLinearNodeDisplacements2D" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::ConstraintLinearNodeDisplacements2D)

void NuTo::ConstraintLinearNodeDisplacements2D::SetNodePtrAfterSerialization(const std::map<uintptr_t, uintptr_t>& mNodeMapCast)
{
    NuTo::ConstraintNode::SetNodePtrAfterSerialization(mNodeMapCast);
}
#endif // ENABLE_SERIALIZATION
