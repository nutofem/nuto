#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION
#include "nuto/mechanics/nodes/NodeGrid1D.h"

//! @brief constructor
NuTo::NodeGrid1D::NodeGrid1D()
{
	mNodeId = 0.;
}

//! @brief constructor
NuTo::NodeGrid1D::NodeGrid1D(int rNodeId)  : NodeBase ()
{
	mNodeId = rNodeId;
}

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void NuTo::NodeGrid1D::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeBase);
    }
#endif // ENABLE_SERIALIZATION


//! @brief returns the number of coordinates of the node
//! @return number of coordinates
int NuTo::NodeGrid1D::GetNumCoordinates()const
{
	return 1;
}

//! @brief set the coordinates
//! @param rCoordinates  given coordinates
void NuTo::NodeGrid1D::SetNodeId(int rNodeId)
{
	mNodeId = rNodeId;
}

//! @brief set the coordinates
//! @param rCoordinates  given coordinates
int NuTo::NodeGrid1D::GetNodeId()const
{
	return mNodeId;
}

//! @brief writes the coordinates of a node to the prescribed pointer
//! @param rCoordinates coordinates
void NuTo::NodeGrid1D::GetCoordinates1D(double rCoordinates[1])const
{
	throw MechanicsException("[NuTo::NodeGrid1D::GetCoordinates1D] to be implemented.");
}

//! @brief sets the global dofs
//! @param rDOF current maximum DOF, this variable is increased within the routine
void NuTo::NodeGrid1D::SetGlobalDofs(int& rDOF)
{
	//empty since coordinates are no DOFs
}

//! @brief write dof values to the node (based on global dof number)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::NodeGrid1D::SetGlobalDofValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
{
	//empty since coordinates are no DOFs
}

//! @brief extract dof values from the node (based on global dof number)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::NodeGrid1D::GetGlobalDofValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
{
	//empty since coordinates are no DOFs
}

//! @brief renumber the global dofs according to predefined ordering
//! @param rMappingInitialToNewOrdering ... mapping from initial ordering to the new ordering
void NuTo::NodeGrid1D::RenumberGlobalDofs(std::vector<int>& rMappingInitialToNewOrdering)
{
	//empty since coordinates are no DOFs
}

//! @brief returns the type of the node
//! @return type
std::string NuTo::NodeGrid1D::GetNodeTypeStr()const
{
	return std::string("NodeGrid1D");
}
