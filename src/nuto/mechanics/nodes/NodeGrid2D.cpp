// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION
#include "nuto/mechanics/nodes/NodeGrid2D.h"

//! @brief constructor
NuTo::NodeGrid2D::NodeGrid2D()
{
	mNodeId = 0.;
}

//! @brief constructor
NuTo::NodeGrid2D::NodeGrid2D(int rNodeId)  : NodeBase ()
{
	mNodeId = rNodeId;
}

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void NuTo::NodeGrid2D::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeBase);
    }
#endif // ENABLE_SERIALIZATION


//! @brief returns the number of coordinates of the node
//! @return number of coordinates
int NuTo::NodeGrid2D::GetNumCoordinates()const
{
	return 2;
}

//! @brief set node id
//! @param rNodeId  given node id
void NuTo::NodeGrid2D::SetNodeId(int rNodeId)
{
	mNodeId = rNodeId;
}

//! @brief get node id
//! @return id given node id
int NuTo::NodeGrid2D::GetNodeId()const
{
	return mNodeId;
}

//! @brief writes the coordinates of a node to the prescribed pointer
//! @param rCoordinates coordinates
void NuTo::NodeGrid2D::GetCoordinates2D(double rCoordinates[2])const
{
	throw MechanicsException("[NuTo::NodeGrid2D::GetCoordinates2D] to be implemented.");
}

//! @brief sets the global dofs
//! @param rDOF current maximum DOF, this variable is increased within the routine
void NuTo::NodeGrid2D::SetGlobalDofs(int& rDOF)
{
	//empty since coordinates are no DOFs
}

//! @brief write dof values to the node (based on global dof number)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::NodeGrid2D::SetGlobalDofValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
{
	//empty since coordinates are no DOFs
}

//! @brief extract dof values from the node (based on global dof number)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::NodeGrid2D::GetGlobalDofValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
{
	//empty since coordinates are no DOFs
}

//! @brief write first time derivative of the dof values (e.g. velocities) to the node (based on global dof number)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::NodeGrid2D::SetGlobalDofFirstTimeDerivativeValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
{
	//empty since coordinates are no DOFs
}

//! @brief extract first time derivative of the dof values (e.g. velocities) from the node (based on global dof number)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::NodeGrid2D::GetGlobalDofFirstTimeDerivativeValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
{
	//empty since coordinates are no DOFs
}

//! @brief write second time derivative of the dof values (e.g. accelerations) to the node (based on global dof number)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::NodeGrid2D::SetGlobalDofSecondTimeDerivativeValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
{
	//empty since coordinates are no DOFs
}

//! @brief extract second time derivative of the dof values (e.g. accelerations) from the node (based on global dof number)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::NodeGrid2D::GetGlobalDofSecondTimeDerivativeValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
{
	//empty since coordinates are no DOFs
}

//! @brief renumber the global dofs according to predefined ordering
//! @param rMappingInitialToNewOrdering ... mapping from initial ordering to the new ordering
void NuTo::NodeGrid2D::RenumberGlobalDofs(std::vector<int>& rMappingInitialToNewOrdering)
{
	//empty since coordinates are no DOFs
}

//! @brief returns the type of the node
//! @return type
std::string NuTo::NodeGrid2D::GetNodeTypeStr()const
{
	return std::string("NodeGrid2D");
}


//! @brief returns the type of node as an enum (all the data stored at the node)
//! @return enum
NuTo::Node::eNodeType NuTo::NodeGrid2D::GetNodeType()const
{
    return Node::NodeGrid2D;
}
