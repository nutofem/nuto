// $Id:$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include "nuto/mechanics/nodes/NodeRadius.h"

//! @brief constructor
NuTo::NodeRadius::NodeRadius() : NodeBase ()
{
	this->mRadius=0;
}

//! @brief constructor
NuTo::NodeRadius::NodeRadius (double rRadius)  : NodeBase ()
{
	this->mRadius=rRadius;
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::NodeRadius::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::NodeRadius::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::NodeRadius::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::NodeRadius::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::NodeRadius::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::NodeRadius::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::NodeRadius::serialize(Archive & ar, const unsigned int version)
{
	ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeBase)
	   & BOOST_SERIALIZATION_NVP(mRadius);
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::NodeRadius)
BOOST_CLASS_TRACKING(NuTo::NodeRadius, track_always)
#endif  // ENABLE_SERIALIZATION

//! @brief returns the number of radii
//! @return number of radii
int NuTo::NodeRadius::GetNumRadius()const
{
	return 1;
}

//! @brief returns the radius of the node
//! @param rRadius ... radius
void NuTo::NodeRadius::GetRadius(double rRadius[1])const
{
	rRadius[0]=mRadius;
}

//! @brief set the radius
//! @param rRadius  given radius
void NuTo::NodeRadius::SetRadius(const double rRadius[1])
{
	mRadius=rRadius[0];
}


//! @brief sets the global dofs
//! @param rDOF current maximum DOF, this variable is increased within the routine
void NuTo::NodeRadius::SetGlobalDofs(int& rDOF)
{
}


//! @brief write dof values to the node (based on global dof number)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::NodeRadius::SetGlobalDofValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
{
}

//! @brief extract dof values from the node (based on global dof number)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::NodeRadius::GetGlobalDofValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
{
}

//! @brief write first time derivative of the dof values (e.g. velocities) to the node (based on global dof number)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::NodeRadius::SetGlobalDofFirstTimeDerivativeValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
{
	throw MechanicsException("[NuTo::NodeRadius::SetGlobalDofFirstTimeDerivativeValues] Node of type " + GetNodeTypeStr() + " has no first time derivative dofs.");
}

//! @brief extract first time derivative of the dof values (e.g. velocities) from the node (based on global dof number)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::NodeRadius::GetGlobalDofFirstTimeDerivativeValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
{
	throw MechanicsException("[NuTo::NodeRadius::GetGlobalDofFirstTimeDerivativeValues] Node of type " + GetNodeTypeStr() + " has no first time derivative dofs.");
}

//! @brief write second time derivative of the dof values (e.g. accelerations) to the node (based on global dof number)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::NodeRadius::SetGlobalDofSecondTimeDerivativeValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
{
	throw MechanicsException("[NuTo::NodeRadius::SetGlobalDofSecondTimeDerivativeValues] Node of type " + GetNodeTypeStr() + " has no second time derivative dofs.");
}

//! @brief extract second time derivative of the dof values (e.g. accelerations) from the node (based on global dof number)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::NodeRadius::GetGlobalDofSecondTimeDerivativeValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
{
	throw MechanicsException("[NuTo::NodeRadius::GetGlobalDofSecondTimeDerivativeValues] Node of type " + GetNodeTypeStr() + " has no second time derivative dofs.");
}

//! @brief renumber the global dofs according to predefined ordering
//! @param rMappingInitialToNewOrdering ... mapping from initial ordering to the new ordering
void NuTo::NodeRadius::RenumberGlobalDofs(std::vector<int>& rMappingInitialToNewOrdering)
{
}

//! @brief returns the type of the node
//! @return type
std::string NuTo::NodeRadius::GetNodeTypeStr()const
{
	return std::string("NodeRadius");
}

//! @brief returns the type of node as an enum (all the data stored at the node)
//! @return enum
NuTo::Node::eNodeType NuTo::NodeRadius::GetNodeType()const
{
    return Node::NodeRadius;
}

