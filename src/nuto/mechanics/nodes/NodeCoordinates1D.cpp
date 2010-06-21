#include <iostream>

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION
#include "nuto/mechanics/nodes/NodeCoordinates1D.h"

//! @brief constructor
NuTo::NodeCoordinates1D::NodeCoordinates1D()
{
	mCoordinates[0] = 0.;
}

//! @brief constructor
NuTo::NodeCoordinates1D::NodeCoordinates1D(const double rCoordinates[1])  : NodeBase ()
{
	mCoordinates[0] = rCoordinates[0];
}

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void NuTo::NodeCoordinates1D::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeBase)
           & BOOST_SERIALIZATION_NVP(mCoordinates);
     }
#endif // ENABLE_SERIALIZATION


//! @brief returns the number of coordinates of the node
//! @return number of coordinates
int NuTo::NodeCoordinates1D::GetNumCoordinates()const
{
	return 1;
}

//! @brief set the coordinates
//! @param rCoordinates  given coordinates
void NuTo::NodeCoordinates1D::SetCoordinates1D(const double rCoordinates[1])
{
	mCoordinates[0] = rCoordinates[0];
}

//! @brief writes the coordinates of a node to the prescribed pointer
//! @param rCoordinates coordinates
void NuTo::NodeCoordinates1D::GetCoordinates1D(double rCoordinates[1])const
{
	rCoordinates[0] = mCoordinates[0];
}

//! @brief returns the coordinate of a given direction of the node
//! @param rIndex index of the direction
//! @return coordinate
double NuTo::NodeCoordinates1D::GetCoordinate(short rIndex)const
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
void NuTo::NodeCoordinates1D::SetGlobalDofs(int& rDOF)
{
	//empty since coordinates are no DOFs
}

//! @brief write dof values to the node (based on global dof number)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::NodeCoordinates1D::SetGlobalDofValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
{
	//empty since coordinates are no DOFs
}

//! @brief extract dof values from the node (based on global dof number)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::NodeCoordinates1D::GetGlobalDofValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
{
	//empty since coordinates are no DOFs
}

//! @brief renumber the global dofs according to predefined ordering
//! @param rMappingInitialToNewOrdering ... mapping from initial ordering to the new ordering
void NuTo::NodeCoordinates1D::RenumberGlobalDofs(std::vector<int>& rMappingInitialToNewOrdering)
{
	//empty since coordinates are no DOFs
}

//! @brief returns the type of the node
//! @return type
std::string NuTo::NodeCoordinates1D::GetNodeTypeStr()const
{
	return std::string("NodeCoordinates1D");
}
