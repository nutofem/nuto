#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION
#include "nuto/mechanics/nodes/NodeBase.h"

//! @brief constructor
NuTo::NodeBase::NodeBase()
{}

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void NuTo::NodeBase::serialize(Archive & ar, const unsigned int version)
    {}
#endif // ENABLE_SERIALIZATION

//! @brief returns the number of coordinates of the node
//! @return number of coordinates
int NuTo::NodeBase::GetNumCoordinates()const
{
	return 0;
}

//! @brief returns the coordinates of the node
//! @return coordinates
void NuTo::NodeBase::GetCoordinates1D(double rCoordinates[1])const
{
	throw MechanicsException("[NuTo::NodeBase::GetCoordinates1D] Node of type " + GetNodeTypeStr() + " has no 1D coordinates.");
}

//! @brief set the coordinates
//! @param rCoordinates  given coordinates
void NuTo::NodeBase::SetCoordinates1D(const double rCoordinates[1])
{
	throw MechanicsException("[NuTo::NodeBase::SetCoordinates1D] Node of type " + GetNodeTypeStr() + " has no 1D coordinates.");
}

//! @brief returns the coordinates of the node
//! @return coordinates
void NuTo::NodeBase::GetCoordinates2D(double rCoordinates[2])const
{
	throw MechanicsException("[NuTo::NodeBase::GetCoordinates1D] Node of type " + GetNodeTypeStr() + " has no 2D coordinates.");
}

//! @brief set the coordinates
//! @param rCoordinates  given coordinates
void NuTo::NodeBase::SetCoordinates2D(const double rCoordinates[1])
{
	throw MechanicsException("[NuTo::NodeBase::SetCoordinates2D] Node of type " + GetNodeTypeStr() + " has no 2D coordinates.");
}

//! @brief returns the coordinates of the node
//! @return coordinates
void NuTo::NodeBase::GetCoordinates3D(double rCoordinates[3])const
{
	throw MechanicsException("[NuTo::NodeBase::GetCoordinates1D] Node of type " + GetNodeTypeStr() + " has no 3D coordinates.");
}

//! @brief set the coordinates
//! @param rCoordinates  given coordinates
void NuTo::NodeBase::SetCoordinates3D(const double rCoordinates[1])
{
	throw MechanicsException("[NuTo::NodeBase::SetCoordinates3D] Node of type " + GetNodeTypeStr() + " has no 3D coordinates.");
}

//! @brief returns the number of coordinates of the node
//! @param rIndex index of the component (e.g. 0,1,2)
//! @return coordinate
double NuTo::NodeBase::GetCoordinate(short rIndex)const
{
	throw MechanicsException("[NuTo::NodeBase::GetCoordinates] Node of type " + GetNodeTypeStr() + " has no coordinates.");
}

//! @brief returns the number of displacements of the node
//! @return number of displacements
int NuTo::NodeBase::GetNumDisplacements()const
{
	return 0;
}

//! @brief gives the global DOF of a displacement component
//! @param rComponent component
//! @return global DOF
int NuTo::NodeBase::GetDofDisplacement(int rComponent)const
{
	throw MechanicsException("[NuTo::NodeBase::GetDofDisplacement] Node of type " + GetNodeTypeStr() + " has no displacements.");
}

//! @brief returns the displacements of the node
//! @return displacement
void NuTo::NodeBase::GetDisplacements1D(double rCoordinates[1])const
{
	throw MechanicsException("[NuTo::NodeBase::GetDisplacements1D] Node of type " + GetNodeTypeStr() + " has no 1D displacements.");
}

//! @brief set the displacements
//! @param rDisplacements  given displacements
void NuTo::NodeBase::SetDisplacements1D(const double rDisplacements[1])
{
	throw MechanicsException("[NuTo::NodeBase::SetDisplacements1D] Node of type " + GetNodeTypeStr() + " has no 1D displacements.");
}

//! @brief returns the displacements of the node
//! @return displacement
void NuTo::NodeBase::GetDisplacements2D(double rCoordinates[2])const
{
	throw MechanicsException("[NuTo::NodeBase::GetDisplacements2D] Node of type " + GetNodeTypeStr() + " has no 2D displacements.");
}

//! @brief set the displacements
//! @param rDisplacements  given displacements
void NuTo::NodeBase::SetDisplacements2D(const double rDisplacements[2])
{
	throw MechanicsException("[NuTo::NodeBase::SetDisplacements2D] Node of type " + GetNodeTypeStr() + " has no 2D displacements.");
}

//! @brief returns the displacements of the node
//! @return displacement
void NuTo::NodeBase::GetDisplacements3D(double rCoordinates[3])const
{
	throw MechanicsException("[NuTo::NodeBase::GetDisplacements3D] Node of type " + GetNodeTypeStr() + " has no 3D displacements.");
}

//! @brief set the displacements
//! @param rDisplacements  given displacements
void NuTo::NodeBase::SetDisplacements3D(const double rDisplacements[3])
{
	throw MechanicsException("[NuTo::NodeBase::SetDisplacements3D] Node of type " + GetNodeTypeStr() + " has no 3D displacements.");
}

//! @brief returns the displacements of the node
//! @return displacement
double NuTo::NodeBase::GetDisplacement(short rIndex)const
{
	throw MechanicsException("[NuTo::NodeBase::GetDisplacement] Node of type " + GetNodeTypeStr() + " has no displacements.");
}

//! @brief returns the number of Rotations of the node
//! @return number of Rotations
int NuTo::NodeBase::GetNumRotations()const
{
	return 0;
}

//! @brief gives the global DOF of a Rotation component
//! @param rComponent component
//! @return global DOF
int NuTo::NodeBase::GetDofRotation(int rComponent)const
{
	throw MechanicsException("[NuTo::NodeBase::GetDofRotation] Node of type " + GetNodeTypeStr() + " has no rotations.");
}

//! @brief returns the Rotations of the node
//! @return Rotation
void NuTo::NodeBase::GetRotations2D(double rCoordinates[2])const
{
	throw MechanicsException("[NuTo::NodeBase::GetRotations2D] Node of type " + GetNodeTypeStr() + " has no 2D rotations.");
}

//! @brief set the Rotations
//! @param rRotations  given Rotations
void NuTo::NodeBase::SetRotations2D(const double rRotations[2])
{
	throw MechanicsException("[NuTo::NodeBase::SetRotations2D] Node of type " + GetNodeTypeStr() + " has no 2D rotations.");
}

//! @brief returns the Rotations of the node
//! @return Rotation
void NuTo::NodeBase::GetRotations3D(double rCoordinates[3])const
{
	throw MechanicsException("[NuTo::NodeBase::GetRotations3D] Node of type " + GetNodeTypeStr() + " has no 3D rotations.");
}

//! @brief set the Rotations
//! @param rRotations  given Rotations
void NuTo::NodeBase::SetRotations3D(const double rRotations[3])
{
	throw MechanicsException("[NuTo::NodeBase::SetRotations3D] Node of type " + GetNodeTypeStr() + " has no 3D rotations.");
}

//! @brief returns the Rotations of the node
//! @return Rotation
double NuTo::NodeBase::GetRotation(short rIndex)const
{
	throw MechanicsException("[NuTo::NodeBase::GetRotation] Node of type " + GetNodeTypeStr() + " has no rotations.");
}

//! @brief returns the number of temperatures of the node
//! @return number of temperatures
int NuTo::NodeBase::GetNumTemperatures()const
{
	return 0;
}

//! @brief gives the global DOF of a temperature component
//! @param rComponent component
//! @return global DOF
int NuTo::NodeBase::GetDofTemperature(int rComponent)const
{
	throw MechanicsException("[NuTo::NodeBase::GetDofTemperature] Node of type " + GetNodeTypeStr() + " has no temperatures.");
}

//! @brief gives the grid number of the node, this is only implemented for the Grid node, since it stores the grid number
//! @return NodeGridNum
int NuTo::NodeBase::GetNodeGridNum()const
{
	throw MechanicsException("[NuTo::NodeBase::GetNodeGridNum] Node of type " + GetNodeTypeStr() + " does not store his id - access the id via the structure.");
}
