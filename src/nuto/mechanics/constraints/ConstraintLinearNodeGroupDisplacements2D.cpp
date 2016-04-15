// $Id$
#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include <iostream>

#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/groups/Group.h"
#include "nuto/mechanics/constraints/ConstraintLinearNodeGroupDisplacements2D.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"

//! @brief constructor
NuTo::ConstraintLinearNodeGroupDisplacements2D::ConstraintLinearNodeGroupDisplacements2D(const Group<NodeBase>* rGroup, const NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rDirection, double rValue) :
        ConstraintNodeGroup(rGroup), ConstraintLinear()
{
    if (rDirection.GetNumColumns()!=1 || rDirection.GetNumRows()!=2)
        throw MechanicsException("[NuTo::ConstraintLinearNodeGroupDisplacements2D::ConstraintLinearNodeGroupDisplacements2D] Dimension of the direction matrix must be equal to the dimension of the structure.");

    memcpy(mDirection,rDirection.data(),2*sizeof(double));
    //normalize the direction
    double norm = sqrt(mDirection[0]*mDirection[0]+mDirection[1]*mDirection[1]);
    if (norm < 1e-14)
    {
        throw MechanicsException("[NuTo::ConstraintLinearNodeGroupDisplacements2D::ConstraintLinearNodeGroupDisplacements2D] direction vector has zero length.");
    }
    double invNorm = 1./norm;
    mDirection[0]*=invNorm;
    mDirection[1]*=invNorm;
    mRHS = rValue;
}

//! @brief returns the number of constraint equations
//! @return number of constraints
int NuTo::ConstraintLinearNodeGroupDisplacements2D::GetNumLinearConstraints()const
{
    return mGroup->GetNumMembers();
}

//!@brief sets/modifies the right hand side of the constraint equation
//!@param rRHS new right hand side
void NuTo::ConstraintLinearNodeGroupDisplacements2D::SetRHS(double rRHS)
{
	mRHS=rRHS;
}

//!@brief returns the right hand side of the constraint equation
double NuTo::ConstraintLinearNodeGroupDisplacements2D::GetRHS()const
{
	return mRHS;
}

//! @brief adds the constraint equations to the matrix
//! @param curConstraintEquation (is incremented during the function call)
//! @param rConstraintMatrix (the first row where a constraint equation is added is given by curConstraintEquation)
void NuTo::ConstraintLinearNodeGroupDisplacements2D::AddToConstraintMatrix(int& curConstraintEquation,
        NuTo::SparseMatrix<double>& rConstraintMatrix)const
{
    for (Group<NodeBase>::const_iterator itNode=mGroup->begin(); itNode!=mGroup->end(); itNode++)
    {
        if (itNode->second->GetNumDisplacements()==0)
            throw MechanicsException(__PRETTY_FUNCTION__,"Node does not have displacements");

        if (std::abs(mDirection[0])>1e-18)
            rConstraintMatrix.AddValue(curConstraintEquation,itNode->second->GetDofDisplacement(0),mDirection[0]);
        if (std::abs(mDirection[1])>1e-18)
            rConstraintMatrix.AddValue(curConstraintEquation,itNode->second->GetDofDisplacement(1),mDirection[1]);

        curConstraintEquation++;
    }
}

//!@brief writes for the current constraint equation(s) the rhs into the vector
// (in case of more than one equation per constraint, curConstraintEquation is increased based on the number of constraint equations per constraint)
//! @param curConstraintEquation (is incremented during the function call)
//! @param rConstraintMatrix (the first row where a constraint equation is added is given by curConstraintEquation)
void NuTo::ConstraintLinearNodeGroupDisplacements2D::GetRHS(int& curConstraintEquation,NuTo::FullVector<double,Eigen::Dynamic>& rRHS)const
{
    for (Group<NodeBase>::const_iterator itNode=mGroup->begin(); itNode!=mGroup->end(); itNode++)
    {
        rRHS(curConstraintEquation,0) = mRHS;
        if (itNode->second->GetNumDisplacements()==0)
            throw MechanicsException("[NuTo::ConstraintLinearNodeGroupDisplacements2D::AddToConstraintMatrix] Node does not have displacements or has more than two displacement components.");

        curConstraintEquation++;
    }
}

#ifdef ENABLE_SERIALIZATION
// serialize
template void NuTo::ConstraintLinearNodeGroupDisplacements2D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinearNodeGroupDisplacements2D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinearNodeGroupDisplacements2D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinearNodeGroupDisplacements2D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinearNodeGroupDisplacements2D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinearNodeGroupDisplacements2D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ConstraintLinearNodeGroupDisplacements2D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ConstraintLinearNodeGroupDisplacements2D" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstraintNodeGroup)
       & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstraintLinear)
       & BOOST_SERIALIZATION_NVP(mRHS)
       & BOOST_SERIALIZATION_NVP(mDirection);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ConstraintLinearNodeGroupDisplacements2D" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::ConstraintLinearNodeGroupDisplacements2D)
#endif // ENABLE_SERIALIZATION
