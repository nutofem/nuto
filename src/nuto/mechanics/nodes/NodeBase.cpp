// $Id$

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
{
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::NodeBase::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::NodeBase::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::NodeBase::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::NodeBase::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::NodeBase::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::NodeBase::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::NodeBase::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize NodeBase" << std::endl;
#endif
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize NodeBase" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::NodeBase)
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::NodeBase)
BOOST_CLASS_TRACKING(NuTo::NodeBase, track_always)
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

//! @brief extract all dof numbers from the node (based on global dof number)
int* NuTo::NodeBase::GetGlobalDofs()
{
	throw MechanicsException("[NuTo::NodeBase::GetGlobalDofs] Node of type " + GetNodeTypeStr() + " dofs.");
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
	throw MechanicsException("[NuTo::NodeBase::SetDofDisplacements3D] Node of type " + GetNodeTypeStr() + " has no 3D displacements.");
}

//! @brief returns the displacements of the node
//! @return displacement
double NuTo::NodeBase::GetDisplacement(short rIndex)const
{
	throw MechanicsException("[NuTo::NodeBase::GetDisplacement] Node of type " + GetNodeTypeStr() + " has no displacements.");
}

//! @brief returns the number of displacements of the node
//! @return number of displacements
int NuTo::NodeBase::GetNumFineScaleDisplacements()const
{
    return 0;
}

//! @brief gives the global DOF of a displacement component
//! @param rComponent component
//! @return global DOF
int NuTo::NodeBase::GetDofFineScaleDisplacement(int rComponent)const
{
    throw MechanicsException("[NuTo::NodeBase::GetDofFineScaleDisplacement] Node of type " + GetNodeTypeStr() + " has no displacements.");
}

//! @brief set the displacements
//! @param rDisplacements  given displacements
void NuTo::NodeBase::SetFineScaleDisplacements2D(const double rDisplacements[2])
{
    throw MechanicsException("[NuTo::NodeBase::SetFineScaleDisplacements2D] Node of type " + GetNodeTypeStr() + " has no fine scale displacements.");
}

//! @brief writes the displacements of a node to the prescribed pointer
//! @param rDisplacements displacements
void NuTo::NodeBase::GetFineScaleDisplacements2D(double rDisplacements[2])const
{
    throw MechanicsException("[NuTo::NodeBase::GetFineScaleDisplacements2D] Node of type " + GetNodeTypeStr() + " has no fine scale displacements.");
}

// returns the number of velocities of the node
int NuTo::NodeBase::GetNumVelocities()const
{
	return 0;
}

// gives the global DOF of a velocity component (Note: the velocity dof number is derived from the displacement dof)
int NuTo::NodeBase::GetDofVelocity(int rComponent)const
{
	throw MechanicsException("[NuTo::NodeBase::GetDofVelocity] Node of type " + GetNodeTypeStr() + " has no velocities.");
}

// returns the velocities of the node
void NuTo::NodeBase::GetVelocities1D(double rVelocities[1]) const
{
	throw MechanicsException("[NuTo::NodeBase::GetVelocities1D] Node of type " + GetNodeTypeStr() + " has no 1D velocities.");
}

//! set the velocities
void NuTo::NodeBase::SetVelocities1D(const double rVelocities[1])
{
	throw MechanicsException("[NuTo::NodeBase::SetVelocities1D] Node of type " + GetNodeTypeStr() + " has no 1D velocities.");
}

// returns the velocities of the node
void NuTo::NodeBase::GetVelocities2D(double rVelocities[2]) const
{
	throw MechanicsException("[NuTo::NodeBase::GetVelocities2D] Node of type " + GetNodeTypeStr() + " has no 2D velocities.");
}

//! set the velocities
void NuTo::NodeBase::SetVelocities2D(const double rVelocities[2])
{
	throw MechanicsException("[NuTo::NodeBase::SetVelocities2D] Node of type " + GetNodeTypeStr() + " has no 2D velocities.");
}

// returns the velocities of the node
void NuTo::NodeBase::GetVelocities3D(double rVelocities[3]) const
{
	throw MechanicsException("[NuTo::NodeBase::GetVelocities3D] Node of type " + GetNodeTypeStr() + " has no 3D velocities.");
}

//! set the velocities
void NuTo::NodeBase::SetVelocities3D(const double rVelocities[3])
{
	throw MechanicsException("[NuTo::NodeBase::SetVelocities3D] Node of type " + GetNodeTypeStr() + " has no 3D velocities.");
}

// returns the number of accelerations of the node
int NuTo::NodeBase::GetNumAccelerations()const
{
	return 0;
}

// gives the global DOF of a acceleration component (Note: the acceleration dof number is derived from the displacement dof)
int NuTo::NodeBase::GetDofAcceleration(int rComponent)const
{
	throw MechanicsException("[NuTo::NodeBase::GetDofAcceleration] Node of type " + GetNodeTypeStr() + " has no accelerations.");
}

// returns the velocities of the node
void NuTo::NodeBase::GetAccelerations1D(double rAccelerations[1]) const
{
	throw MechanicsException("[NuTo::NodeBase::GetAccelerations1D] Node of type " + GetNodeTypeStr() + " has no 1D accelerations.");
}

//! set the velocities
void NuTo::NodeBase::SetAccelerations1D(const double rAccelerations[1])
{
	throw MechanicsException("[NuTo::NodeBase::SetAccelerations1D] Node of type " + GetNodeTypeStr() + " has no 1D accelerations.");
}

// returns the velocities of the node
void NuTo::NodeBase::GetAccelerations2D(double rAccelerations[2]) const
{
	throw MechanicsException("[NuTo::NodeBase::GetAccelerations2D] Node of type " + GetNodeTypeStr() + " has no 2D accelerations.");
}

//! set the velocities
void NuTo::NodeBase::SetAccelerations2D(const double rAccelerations[2])
{
	throw MechanicsException("[NuTo::NodeBase::SetAccelerations2D] Node of type " + GetNodeTypeStr() + " has no 2D accelerations.");
}

// returns the velocities of the node
void NuTo::NodeBase::GetAccelerations3D(double rAccelerations[3]) const
{
	throw MechanicsException("[NuTo::NodeBase::GetAccelerations3D] Node of type " + GetNodeTypeStr() + " has no 3D accelerations.");
}

//! set the velocities
void NuTo::NodeBase::SetAccelerations3D(const double rAccelerations[3])
{
	throw MechanicsException("[NuTo::NodeBase::SetAccelerations3D] Node of type " + GetNodeTypeStr() + " has no 3D accelerations.");
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
