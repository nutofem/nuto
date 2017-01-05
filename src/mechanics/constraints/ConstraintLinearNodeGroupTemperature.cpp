// $Id: ConstraintLinearNodeGroupTemperature.cpp 596 2012-03-02 21:11:48Z unger3 $

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include <iostream>

#include "mechanics/MechanicsException.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/groups/Group.h"
#include "mechanics/constraints/ConstraintLinearNodeGroupTemperature.h"
#include "math/FullMatrix.h"
#include "math/SparseMatrixCSRGeneral.h"

//! @brief constructor
NuTo::ConstraintLinearNodeGroupTemperature::ConstraintLinearNodeGroupTemperature(const Group<NodeBase>* rGroup, double rValue) :
        ConstraintNodeGroup(rGroup), ConstraintLinear()
{
    mRHS = rValue;
}

//! @brief returns the number of constraint equations
//! @return number of constraints
int NuTo::ConstraintLinearNodeGroupTemperature::GetNumLinearConstraints()const
{
    return mGroup->GetNumMembers();
}

//!@brief sets/modifies the right hand side of the constraint equation
//!@param rRHS new right hand side
void NuTo::ConstraintLinearNodeGroupTemperature::SetRHS(double rRHS)
{
	mRHS=rRHS;
}

//!@brief returns the right hand side of the constraint equation
double NuTo::ConstraintLinearNodeGroupTemperature::GetRHS()const
{
	return mRHS;
}

//! @brief adds the constraint equations to the matrix
//! @param curConstraintEquation (is incremented during the function call)
//! @param rConstraintMatrix (the first row where a constraint equation is added is given by curConstraintEquation)
void NuTo::ConstraintLinearNodeGroupTemperature::AddToConstraintMatrix(int& curConstraintEquation,
        NuTo::SparseMatrix<double>& rConstraintMatrix)const
{
    for (Group<NodeBase>::const_iterator itNode=mGroup->begin(); itNode!=mGroup->end(); itNode++)
    {
        if (itNode->second->GetNum(Node::eDof::TEMPERATURE)!=1)
            throw MechanicsException("[NuTo::ConstraintLinearNodeGroupTemperature::AddToConstraintMatrix] Node does not have temperature dof.");

        rConstraintMatrix.AddValue(curConstraintEquation,itNode->second->GetDof(Node::eDof::TEMPERATURE),1);

        curConstraintEquation++;
    }
}

//!@brief writes for the current constraint equation(s) the rhs into the vector
// (in case of more than one equation per constraint, curConstraintEquation is increased based on the number of constraint equations per constraint)
//! @param curConstraintEquation (is incremented during the function call)
//! @param rConstraintMatrix (the first row where a constraint equation is added is given by curConstraintEquation)
void NuTo::ConstraintLinearNodeGroupTemperature::GetRHS(int& curConstraintEquation,Eigen::VectorXd& rRHS)const
{
	for (Group<NodeBase>::const_iterator itNode=mGroup->begin(); itNode!=mGroup->end(); itNode++)
	{
        if (itNode->second->GetNum(Node::eDof::TEMPERATURE)!=1)
			throw MechanicsException("[NuTo::ConstraintLinearNodeGroupTemperature::GetRHS] Node does not have a temperature component.");
		rRHS(curConstraintEquation) = mRHS;
		curConstraintEquation++;
	}
}

NuTo::Node::eDof NuTo::ConstraintLinearNodeGroupTemperature::GetDofType() const
{
    return Node::eDof::TEMPERATURE;
}


#ifdef ENABLE_SERIALIZATION
// serialize
template void NuTo::ConstraintLinearNodeGroupTemperature::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinearNodeGroupTemperature::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinearNodeGroupTemperature::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinearNodeGroupTemperature::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinearNodeGroupTemperature::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinearNodeGroupTemperature::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ConstraintLinearNodeGroupTemperature::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ConstraintLinearNodeGroupTemperature" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstraintNodeGroup)
       & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstraintLinear)
       & BOOST_SERIALIZATION_NVP(mRHS);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ConstraintLinearNodeGroupTemperature" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::ConstraintLinearNodeGroupTemperature)
#endif // ENABLE_SERIALIZATION
