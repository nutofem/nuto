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

//! @brief sets the global dofs
//! @param rDOF current maximum DOF, this variable is increased within the routine
void NuTo::NodeBase::SetGlobalDofs(int& rDOF)
{
	throw MechanicsException("[NuTo::NodeBase::SetGlobalDofs] Node of type " + GetNodeTypeStr() + " dofs.");
}

//! @brief write dof values to the node (based on global dof number)
//! @param rTimeDerivative ... time derivative (e.g. 0 disp, 1 vel, 2 acc)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::NodeBase::SetGlobalDofValues(int rTimeDerivative, const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
{
	throw MechanicsException("[NuTo::NodeBase::GetGlobalDofValues] Node of type " + GetNodeTypeStr() + " dofs.");
}

//! @brief extract dof values from the node (based on global dof number)
//! @param rTimeDerivative ... time derivative (e.g. 0 disp, 1 vel, 2 acc)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::NodeBase::GetGlobalDofValues(int rTimeDerivative, FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
{
	throw MechanicsException("[NuTo::NodeBase::GetGlobalDofValues] Node of type " + GetNodeTypeStr() + " dofs.");
}

//! @brief extract all dof numbers from the node (based on global dof number)
//int* NuTo::NodeBase::GetGlobalDofs()
//{
//	throw MechanicsException("[NuTo::NodeBase::GetGlobalDofs] Node of type " + GetNodeTypeStr() + " dofs.");
//}

//! @brief renumber the global dofs according to predefined ordering
//! @param rMappingInitialToNewOrdering ... mapping from initial ordering to the new ordering
void NuTo::NodeBase::RenumberGlobalDofs(std::vector<int>& rMappingInitialToNewOrdering)
{
	throw MechanicsException("[NuTo::NodeBase::RenumberGlobalDofs] Node of type " + GetNodeTypeStr() + " dofs.");
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

//! @brief returns the number of time derivatives stored at the node
//! @return number of derivatives
int NuTo::NodeBase::GetNumTimeDerivatives()const
{
	throw MechanicsException("[NuTo::NodeBase::GetNumTimeDerivatives] Node of type " + GetNodeTypeStr() + " has no dofs");
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

//! @brief returns the displacements of the node
//! @param rTimeDerivative time derivative (0
//! @return displacement
void NuTo::NodeBase::GetDisplacements1D(int rTimeDerivative, double rCoordinates[1])const
{
	throw MechanicsException("[NuTo::NodeBase::GetDisplacements1D] Node of type " + GetNodeTypeStr() + " has no 1D displacements.");
}

//! @brief set the displacements
//! @param rDisplacements  given displacements
void NuTo::NodeBase::SetDisplacements1D(const double rDisplacements[1])
{
	throw MechanicsException("[NuTo::NodeBase::SetDisplacements1D] Node of type " + GetNodeTypeStr() + " has no 1D displacements.");
}

//! @brief set the displacements
//! @param rTimeDerivative time derivative (0
//! @param rDisplacements  given displacements
void NuTo::NodeBase::SetDisplacements1D(int rTimeDerivative, const double rDisplacements[1])
{
	throw MechanicsException("[NuTo::NodeBase::SetDisplacements1D] Node of type " + GetNodeTypeStr() + " has no 1D displacements.");
}

//! @brief returns the displacements of the node
//! @return displacement
void NuTo::NodeBase::GetDisplacements2D(double rCoordinates[2])const
{
	throw MechanicsException("[NuTo::NodeBase::GetDisplacements2D] Node of type " + GetNodeTypeStr() + " has no 2D displacements.");
}

//! @brief returns the displacements of the node
//! @param rTimeDerivative time derivative (0
//! @return displacement
void NuTo::NodeBase::GetDisplacements2D(int rTimeDerivative, double rCoordinates[2])const
{
	throw MechanicsException("[NuTo::NodeBase::GetDisplacements2D] Node of type " + GetNodeTypeStr() + " has no 2D displacements.");
}

//! @brief set the displacements
//! @param rDisplacements  given displacements
void NuTo::NodeBase::SetDisplacements2D(const double rDisplacements[2])
{
	throw MechanicsException("[NuTo::NodeBase::SetDisplacements2D] Node of type " + GetNodeTypeStr() + " has no 2D displacements.");
}

//! @brief set the displacements
//! @param rTimeDerivative time derivative (0
//! @param rDisplacements  given displacements
void NuTo::NodeBase::SetDisplacements2D(int rTimeDerivative, const double rDisplacements[2])
{
	throw MechanicsException("[NuTo::NodeBase::SetDisplacements2D] Node of type " + GetNodeTypeStr() + " has no 2D displacements.");
}

//! @brief returns the displacements of the node
//! @return displacement
void NuTo::NodeBase::GetDisplacements3D(double rCoordinates[3])const
{
	throw MechanicsException("[NuTo::NodeBase::GetDisplacements3D] Node of type " + GetNodeTypeStr() + " has no 3D displacements.");
}

//! @brief returns the displacements of the node
//! @param rTimeDerivative time derivative (0
//! @return displacement
void NuTo::NodeBase::GetDisplacements3D(int rTimeDerivative, double rCoordinates[3])const
{
	throw MechanicsException("[NuTo::NodeBase::GetDisplacements3D] Node of type " + GetNodeTypeStr() + " has no 3D displacements.");
}

//! @brief set the displacements
//! @param rDisplacements  given displacements
void NuTo::NodeBase::SetDisplacements3D(const double rDisplacements[3])
{
	throw MechanicsException("[NuTo::NodeBase::SetDofDisplacements3D] Node of type " + GetNodeTypeStr() + " has no 3D displacements.");
}

//! @brief set the displacements
//! @param rTimeDerivative time derivative (0
//! @param rDisplacements  given displacements
void NuTo::NodeBase::SetDisplacements3D(int rTimeDerivative, const double rDisplacements[3])
{
	throw MechanicsException("[NuTo::NodeBase::SetDofDisplacements3D] Node of type " + GetNodeTypeStr() + " has no 3D displacements.");
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

//! @brief returns the Rotations of the node
//! @param rTimeDerivative time derivative (0
//! @return Rotation
void NuTo::NodeBase::GetRotations2D(int rTimeDerivative, double rCoordinates[2])const
{
	throw MechanicsException("[NuTo::NodeBase::GetRotations2D] Node of type " + GetNodeTypeStr() + " has no 2D rotations.");
}

//! @brief set the Rotations
//! @param rRotations  given Rotations
void NuTo::NodeBase::SetRotations2D(const double rRotations[2])
{
	throw MechanicsException("[NuTo::NodeBase::SetRotations2D] Node of type " + GetNodeTypeStr() + " has no 2D rotations.");
}

//! @brief set the Rotations
//! @param rTimeDerivative time derivative (0
//! @param rRotations  given Rotations
void NuTo::NodeBase::SetRotations2D(int rTimeDerivative, const double rRotations[2])
{
	throw MechanicsException("[NuTo::NodeBase::SetRotations2D] Node of type " + GetNodeTypeStr() + " has no 2D rotations.");
}

//! @brief returns the Rotations of the node
//! @return Rotation
void NuTo::NodeBase::GetRotations3D(double rCoordinates[3])const
{
	throw MechanicsException("[NuTo::NodeBase::GetRotations3D] Node of type " + GetNodeTypeStr() + " has no 3D rotations.");
}

//! @brief returns the Rotations of the node
//! @param rTimeDerivative time derivative (0
//! @return Rotation
void NuTo::NodeBase::GetRotations3D(int rTimeDerivative, double rCoordinates[3])const
{
	throw MechanicsException("[NuTo::NodeBase::GetRotations3D] Node of type " + GetNodeTypeStr() + " has no 3D rotations.");
}

//! @brief set the Rotations
//! @param rRotations  given Rotations
void NuTo::NodeBase::SetRotations3D(const double rRotations[3])
{
	throw MechanicsException("[NuTo::NodeBase::SetRotations3D] Node of type " + GetNodeTypeStr() + " has no 3D rotations.");
}

//! @brief set the Rotations
//! @param rTimeDerivative time derivative (0
//! @param rRotations  given Rotations
void NuTo::NodeBase::SetRotations3D(int rTimeDerivative, const double rRotations[3])
{
	throw MechanicsException("[NuTo::NodeBase::SetRotations3D] Node of type " + GetNodeTypeStr() + " has no 3D rotations.");
}

//! @brief returns the Rotations of the node
//! @return Rotation
double NuTo::NodeBase::GetRotation(short rIndex)const
{
	throw MechanicsException("[NuTo::NodeBase::GetRotation] Node of type " + GetNodeTypeStr() + " has no rotations.");
}


//! @brief returns the number of radii
//! @return number of radii
int NuTo::NodeBase::GetNumRadius()const
{
	return 0;
}

//! @brief returns the radius of the node
//! @param rRadius ... radius
void NuTo::NodeBase::GetRadius(double rRadius[1])const
{
	throw MechanicsException("[NuTo::NodeBase::GetRadius] Node of type " + GetNodeTypeStr() + " has no radius.");
}

//! @brief set the radius
//! @param rRadius  given radius
void NuTo::NodeBase::SetRadius(const double rRadius[1])
{
	throw MechanicsException("[NuTo::NodeBase::SetRadius] Node of type " + GetNodeTypeStr() + " has no radius.");
}

//! @brief returns the number of temperatures of the node
//! @return number of temperatures
int NuTo::NodeBase::GetNumTemperatures()const
{
	return 0;
}

//! @brief returns the temperature of the node
//! @return temperature
void NuTo::NodeBase::SetTemperature(const double rTemperature[1])
{
	throw MechanicsException("[NuTo::NodeBase::SetTemperature] Node of type " + GetNodeTypeStr() + " has no temperatures.");
}

//! @brief returns the temperature of the node
//! @return temperature
void NuTo::NodeBase::SetTemperature(int rTimeDerivative, const double rTemperature[1])
{
	throw MechanicsException("[NuTo::NodeBase::SetTemperature] Node of type " + GetNodeTypeStr() + " has no temperatures.");
}

//! @brief returns the temperature of the node
//! @return temperature
void NuTo::NodeBase::GetTemperature(double rTemperature[1])const
{
	throw MechanicsException("[NuTo::NodeBase::GetTemperature] Node of type " + GetNodeTypeStr() + " has no temperatures.");
}

//! @brief returns the temperature of the node
//! @return temperature
void NuTo::NodeBase::GetTemperature(int rTimeDerivative, double rTemperature[1])const
{
	throw MechanicsException("[NuTo::NodeBase::GetTemperature] Node of type " + GetNodeTypeStr() + " has no temperatures.");
}

//! @brief gives the global DOF of a temperature component
//! @param rComponent component
//! @return global DOF
int NuTo::NodeBase::GetDofTemperature()const
{
	throw MechanicsException("[NuTo::NodeBase::GetDofTemperature] Node of type " + GetNodeTypeStr() + " has no temperatures.");
}


#ifdef ENABLE_VISUALIZE
void NuTo::NodeBase::Visualize(VisualizeUnstructuredGrid& rVisualize, const boost::ptr_list<NuTo::VisualizeComponentBase>& rWhat) const
{
	double coordinates[3]={0,0,0};
	switch (this->GetNumCoordinates())
	{
	case 1:
		this->GetCoordinates1D(coordinates);
		break;
	case 2:
		this->GetCoordinates2D(coordinates);
		break;
	case 3:
		this->GetCoordinates3D(coordinates);
		break;
	default:
		throw MechanicsException("[NuTo::NodeBase::Visualize] node has neither coordinates in 1D, 2D or 3D.");
	}
	unsigned int PointId = rVisualize.AddPoint(coordinates);

    // store data
    boost::ptr_list<VisualizeComponentBase>::const_iterator WhatIter = rWhat.begin();
    while (WhatIter != rWhat.end())
    {
        switch (WhatIter->GetComponentEnum())
        {
			case NuTo::VisualizeBase::DISPLACEMENTS:
			{
				double displacements[3]={0,0,0};
				switch (this->GetNumDisplacements())
				{
				case 1:
					this->GetDisplacements1D(displacements);
					break;
				case 2:
					this->GetDisplacements2D(displacements);
					break;
				case 3:
					this->GetDisplacements3D(displacements);
					break;
				default:
					throw MechanicsException("[NuTo::NodeBase::Visualize] node has neither displacements in 1D, 2D or 3D.");
				}
					rVisualize.SetPointDataVector(PointId, WhatIter->GetComponentName(), displacements);
			}
				break;
			case NuTo::VisualizeBase::ROTATION:
			{
				double rotations[3]={0,0,0};
				switch (this->GetNumRotations())
				{
				case 1:
					this->GetRotations2D(rotations);
					break;
				case 3:
					this->GetRotations3D(rotations);
					break;
				default:
					throw MechanicsException("[NuTo::NodeBase::Visualize] node has neither rotations in 2D or 3D.");
				}
					rVisualize.SetPointDataVector(PointId, WhatIter->GetComponentName(), rotations);
			}
				break;
			case NuTo::VisualizeBase::VELOCITY:
			{
				double velocities[3]={0,0,0};
				switch (this->GetNumDisplacements())
				{
				case 1:
					this->GetDisplacements1D(1,velocities);
					break;
				case 2:
					this->GetDisplacements2D(1,velocities);
					break;
				case 3:
					this->GetDisplacements3D(1,velocities);
					break;
				default:
					throw MechanicsException("[NuTo::NodeBase::Visualize] node has neither velocities in 1D, 2D or 3D.");
				}
					rVisualize.SetPointDataVector(PointId, WhatIter->GetComponentName(), velocities);
			}
				break;
			case NuTo::VisualizeBase::ANGULAR_VELOCITY:
			{
				double angularVelocities[3]={0,0,0};
				switch (this->GetNumRotations())
				{
				case 1:
					this->GetRotations2D(0,angularVelocities);
					break;
				case 3:
					this->GetRotations3D(0,angularVelocities);
					break;
				default:
					throw MechanicsException("[NuTo::NodeBase::Visualize] node has neither angular velocities in 2D or 3D.");
				}
					rVisualize.SetPointDataVector(PointId, WhatIter->GetComponentName(), angularVelocities);
			}
				break;
			case NuTo::VisualizeBase::ACCELERATION:
			{
				double accelerations[3]={0,0,0};
				switch (this->GetNumDisplacements())
				{
				case 1:
					this->GetDisplacements1D(2,accelerations);
					break;
				case 2:
					this->GetDisplacements2D(2,accelerations);
					break;
				case 3:
					this->GetDisplacements3D(2,accelerations);
					break;
				default:
					throw MechanicsException("[NuTo::NodeBase::Visualize] node has neither accelerations in 1D, 2D or 3D.");
				}
					rVisualize.SetPointDataVector(PointId, WhatIter->GetComponentName(), accelerations);
			}
				break;
			case NuTo::VisualizeBase::ANGULAR_ACCELERATION:
			{
				double angularAccelerations[3]={0,0,0};
				switch (this->GetNumRotations())
				{
				case 1:
					this->GetRotations2D(2,angularAccelerations);
					break;
				case 3:
					this->GetRotations3D(2,angularAccelerations);
					break;
				default:
					throw MechanicsException("[NuTo::NodeBase::Visualize] node has neither angular accelerations in 2D or 3D.");
				}
					rVisualize.SetPointDataVector(PointId, WhatIter->GetComponentName(), angularAccelerations);
			}
				break;
			case NuTo::VisualizeBase::PARTICLE_RADIUS:
			{
				double radius;
				GetRadius(&radius);
				rVisualize.SetPointDataScalar(PointId, WhatIter->GetComponentName(), radius);
			}
				break;
			default:
				;
        }
        WhatIter++;
    }
}
#endif // ENABLE_VISUALIZE

