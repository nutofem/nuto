#include <sstream>

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION
#include "nuto/mechanics/nodes/NodeCoordinates3D.h"
#include "nuto/math/FullMatrix.h"

//! @brief constructor
NuTo::NodeCoordinates3D::NodeCoordinates3D()
{
	mCoordinates[0] = 0.;
	mCoordinates[1] = 0.;
	mCoordinates[2] = 0.;
}

//! @brief constructor
NuTo::NodeCoordinates3D::NodeCoordinates3D(const double rCoordinates[3])  : NodeBase ()
{
	mCoordinates[0] = rCoordinates[0];
	mCoordinates[1] = rCoordinates[1];
	mCoordinates[2] = rCoordinates[2];
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::NodeCoordinates3D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::NodeCoordinates3D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::NodeCoordinates3D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::NodeCoordinates3D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::NodeCoordinates3D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::NodeCoordinates3D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::NodeCoordinates3D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize NodeCoordinates3D" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeBase)
       & BOOST_SERIALIZATION_NVP(mCoordinates);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize NodeCoordinates3D" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::NodeCoordinates3D)
BOOST_CLASS_TRACKING(NuTo::NodeCoordinates3D, track_always)
#endif // ENABLE_SERIALIZATION


//! @brief returns the number of coordinates of the node
//! @return number of coordinates
int NuTo::NodeCoordinates3D::GetNumCoordinates()const
{
	return 3;
}

//! @brief set the coordinates
//! @param rCoordinates  given coordinates
void NuTo::NodeCoordinates3D::SetCoordinates3D(const double rCoordinates[3])
{
	mCoordinates[0] = rCoordinates[0];
	mCoordinates[1] = rCoordinates[1];
	mCoordinates[2] = rCoordinates[2];
}

//! @brief writes the coordinates of a node to the prescribed pointer
//! @param rCoordinates coordinates
void NuTo::NodeCoordinates3D::GetCoordinates3D(double rCoordinates[3])const
{
	rCoordinates[0] = mCoordinates[0];
	rCoordinates[1] = mCoordinates[1];
	rCoordinates[2] = mCoordinates[2];
}

//! @brief returns the coordinate of a given direction of the node
//! @param rIndex index of the direction
//! @return coordinate
double NuTo::NodeCoordinates3D::GetCoordinate(short rIndex)const
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
void NuTo::NodeCoordinates3D::SetGlobalDofs(int& rDOF)
{
	//empty since coordinates are no DOFs
}

//! @brief write dof values to the node (based on global dof number)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::NodeCoordinates3D::SetGlobalDofValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
{
	//empty since coordinates are no DOFs
}

//! @brief extract dof values from the node (based on global dof number)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::NodeCoordinates3D::GetGlobalDofValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
{
	//empty since coordinates are no DOFs
}

//! @brief write first time derivative of the dof values (e.g. velocities) to the node (based on global dof number)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::NodeCoordinates3D::SetGlobalDofFirstTimeDerivativeValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
{
	//empty since coordinates are no DOFs
}

//! @brief extract first time derivative of the dof values (e.g. velocities) from the node (based on global dof number)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::NodeCoordinates3D::GetGlobalDofFirstTimeDerivativeValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
{
	//empty since coordinates are no DOFs
}

//! @brief write second time derivative of the dof values (e.g. accelerations) to the node (based on global dof number)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::NodeCoordinates3D::SetGlobalDofSecondTimeDerivativeValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
{
	//empty since coordinates are no DOFs
}

//! @brief extract second time derivative of the dof values (e.g. accelerations) from the node (based on global dof number)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::NodeCoordinates3D::GetGlobalDofSecondTimeDerivativeValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
{
	//empty since coordinates are no DOFs
}

//! @brief renumber the global dofs according to predefined ordering
//! @param rMappingInitialToNewOrdering ... mapping from initial ordering to the new ordering
void NuTo::NodeCoordinates3D::RenumberGlobalDofs(std::vector<int>& rMappingInitialToNewOrdering)
{
	//empty since coordinates are no DOFs
}

//! @brief returns the type of the node
//! @return type
std::string NuTo::NodeCoordinates3D::GetNodeTypeStr()const
{
	return std::string("NodeCoordinates3D");
}
