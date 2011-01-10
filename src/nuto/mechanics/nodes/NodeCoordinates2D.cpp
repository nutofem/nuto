// $Id$

#include <sstream>

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION
#include "nuto/mechanics/nodes/NodeCoordinates2D.h"
#include "nuto/math/FullMatrix.h"

//! @brief constructor
NuTo::NodeCoordinates2D::NodeCoordinates2D()
{
	mCoordinates[0] = 0.;
    mCoordinates[1] = 0.;
}

//! @brief constructor
NuTo::NodeCoordinates2D::NodeCoordinates2D(const double rCoordinates[2])  : NodeBase ()
{
	mCoordinates[0] = rCoordinates[0];
	mCoordinates[1] = rCoordinates[1];
}


#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::NodeCoordinates2D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::NodeCoordinates2D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::NodeCoordinates2D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::NodeCoordinates2D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::NodeCoordinates2D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::NodeCoordinates2D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::NodeCoordinates2D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize NodeCoordinates2D" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeBase)
       & BOOST_SERIALIZATION_NVP(mCoordinates);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize NodeCoordinates2D" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::NodeCoordinates2D)
BOOST_CLASS_TRACKING(NuTo::NodeCoordinates2D, track_always)
#endif // ENABLE_SERIALIZATION


//! @brief returns the number of coordinates of the node
//! @return number of coordinates
int NuTo::NodeCoordinates2D::GetNumCoordinates()const
{
	return 2;
}

//! @brief set the coordinates
//! @param rCoordinates  given coordinates
void NuTo::NodeCoordinates2D::SetCoordinates2D(const double rCoordinates[2])
{
	mCoordinates[0] = rCoordinates[0];
	mCoordinates[1] = rCoordinates[1];
}

//! @brief writes the coordinates of a node to the prescribed pointer
//! @param rCoordinates coordinates
void NuTo::NodeCoordinates2D::GetCoordinates2D(double rCoordinates[2])const
{
    rCoordinates[0] = mCoordinates[0];
	rCoordinates[1] = mCoordinates[1];
}

//! @brief returns the coordinate of a given direction of the node
//! @param rIndex index of the direction
//! @return coordinate
double NuTo::NodeCoordinates2D::GetCoordinate(short rIndex)const
{
	if((rIndex > this->GetNumCoordinates()) | (rIndex < 0 ) )
	{
		std::stringstream ss;
		ss << rIndex;
		throw NuTo::MechanicsException(
				"[NuTo::NodeCoordinates::GetCoordinate] Error getting coordinate(" + ss.str() + ") for " + GetNodeTypeStr() + ".");
	}

	return mCoordinates[rIndex];
}

//! @brief sets the global dofs
//! @param rDOF current maximum DOF, this variable is increased within the routine
void NuTo::NodeCoordinates2D::SetGlobalDofs(int& rDOF)
{
	//empty since coordinates are no DOFs
}

//! @brief write dof values to the node (based on global dof number)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::NodeCoordinates2D::SetGlobalDofValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
{
	//empty since coordinates are no DOFs
}

//! @brief extract dof values from the node (based on global dof number)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::NodeCoordinates2D::GetGlobalDofValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
{
	//empty since coordinates are no DOFs
}

//! @brief renumber the global dofs according to predefined ordering
//! @param rMappingInitialToNewOrdering ... mapping from initial ordering to the new ordering
void NuTo::NodeCoordinates2D::RenumberGlobalDofs(std::vector<int>& rMappingInitialToNewOrdering)
{
	//empty since coordinates are no DOFs
}

//! @brief write first time derivative of the dof values (e.g. velocities) to the node (based on global dof number)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::NodeCoordinates2D::SetGlobalDofFirstTimeDerivativeValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
{
	//empty since coordinates are no DOFs
}

//! @brief extract first time derivative of the dof values (e.g. velocities) from the node (based on global dof number)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::NodeCoordinates2D::GetGlobalDofFirstTimeDerivativeValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
{
	//empty since coordinates are no DOFs
}

//! @brief write second time derivative of the dof values (e.g. accelerations) to the node (based on global dof number)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::NodeCoordinates2D::SetGlobalDofSecondTimeDerivativeValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
{
	//empty since coordinates are no DOFs
}

//! @brief extract second time derivative of the dof values (e.g. accelerations) from the node (based on global dof number)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::NodeCoordinates2D::GetGlobalDofSecondTimeDerivativeValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
{
	//empty since coordinates are no DOFs
}

//! @brief returns the type of the node
//! @return type
std::string NuTo::NodeCoordinates2D::GetNodeTypeStr()const
{
	return std::string("NodeCoordinates2D");
}


//! @brief returns the type of node as an enum (all the data stored at the node)
//! @return enum
NuTo::Node::eNodeType NuTo::NodeCoordinates2D::GetNodeType()const
{
    return NuTo::Node::NodeCoordinates2D;
}
