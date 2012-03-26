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

// returns the number of velocities of the node
int NuTo::NodeBase::GetNumAngularVelocities()const
{
	return 0;
}

// gives the global DOF of a velocity component (Note: the velocity dof number is derived from the displacement dof)
int NuTo::NodeBase::GetDofAngularVelocity(int rComponent)const
{
	throw MechanicsException("[NuTo::NodeBase::GetDofAngularVelocity] Node of type " + GetNodeTypeStr() + " has no angular velocities.");
}

// returns the velocities of the node
void NuTo::NodeBase::GetAngularVelocities2D(double rVelocities[2]) const
{
	throw MechanicsException("[NuTo::NodeBase::GetAngularVelocities2D] Node of type " + GetNodeTypeStr() + " has no 2D angular velocities.");
}

//! set the velocities
void NuTo::NodeBase::SetAngularVelocities2D(const double rVelocities[2])
{
	throw MechanicsException("[NuTo::NodeBase::SetAngularVelocities2D] Node of type " + GetNodeTypeStr() + " has no 2D angular velocities.");
}

// returns the velocities of the node
void NuTo::NodeBase::GetAngularVelocities3D(double rVelocities[3]) const
{
	throw MechanicsException("[NuTo::NodeBase::GetAngularVelocities3D] Node of type " + GetNodeTypeStr() + " has no 3D angular velocities.");
}

//! set the velocities
void NuTo::NodeBase::SetAngularVelocities3D(const double rVelocities[3])
{
	throw MechanicsException("[NuTo::NodeBase::SetAngularVelocities3D] Node of type " + GetNodeTypeStr() + " has no 3D angular velocities.");
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

// returns the number of accelerations of the node
int NuTo::NodeBase::GetNumAngularAccelerations()const
{
	return 0;
}

// gives the global DOF of a acceleration component (Note: the angular acceleration dof number is derived from the rotation dof)
int NuTo::NodeBase::GetDofAngularAcceleration(int rComponent)const
{
	throw MechanicsException("[NuTo::NodeBase::GetDofAngularAcceleration] Node of type " + GetNodeTypeStr() + " has no angular accelerations.");
}

// returns the angular accelerations of the node
void NuTo::NodeBase::GetAngularAccelerations2D(double rAngularAccelerations[2]) const
{
	throw MechanicsException("[NuTo::NodeBase::GetAngularAccelerations2D] Node of type " + GetNodeTypeStr() + " has no 2D angular accelerations.");
}

//! set the angular accelerations
void NuTo::NodeBase::SetAngularAccelerations2D(const double rAngularAccelerations[2])
{
	throw MechanicsException("[NuTo::NodeBase::SetAngularAccelerations2D] Node of type " + GetNodeTypeStr() + " has no 2D angular accelerations.");
}

// returns the angular accelerations of the node
void NuTo::NodeBase::GetAngularAccelerations3D(double rAngularAccelerations[3]) const
{
	throw MechanicsException("[NuTo::NodeBase::GetAngularAccelerations3D] Node of type " + GetNodeTypeStr() + " has no 3D angular accelerations.");
}

//! set the angular accelerations
void NuTo::NodeBase::SetAngularAccelerations3D(const double rAngularAccelerations[3])
{
	throw MechanicsException("[NuTo::NodeBase::SetAngularAccelerations3D] Node of type " + GetNodeTypeStr() + " has no 3D angular accelerations.");
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

//! @brief gives the global DOF of a temperature component
//! @param rComponent component
//! @return global DOF
int NuTo::NodeBase::GetDofTemperature(int rComponent)const
{
	throw MechanicsException("[NuTo::NodeBase::GetDofTemperature] Node of type " + GetNodeTypeStr() + " has no temperatures.");
}

//! @brief set the shape functions based on the actual oscillations
//! @parameter shape function number (0..2)
void NuTo::NodeBase::SetShapeFunctionMultiscalePeriodic(int rShapeFunction)
{
	throw MechanicsException("[NuTo::NodeBase::SetShapeFunctionMultiscalePeriodic] Node of type " + GetNodeTypeStr() + " is not a multiscale node.");
}

//! @brief returns the shape function for the periodic bc for the nodes
const boost::array<double,3>& NuTo::NodeBase::GetShapeFunctionMultiscalePeriodicX()const
{
	throw MechanicsException("[NuTo::NodeBase::GetShapeFunctionMultiscalePeriodicX] Node of type " + GetNodeTypeStr() + " is not a multiscale node.");
}

//! @brief returns the shape function for the periodic bc for the nodes
const boost::array<double,3>& NuTo::NodeBase::GetShapeFunctionMultiscalePeriodicY() const
{
	throw MechanicsException("[NuTo::NodeBase::GetShapeFunctionMultiscalePeriodicY] Node of type " + GetNodeTypeStr() + " is not a multiscale node.");
}

//! @brief scales the shape functions
//! @parameter rShapeFunction  (1..3 corresponding to macro strains exx, eyy, and gxy)
//! @parameter rScalingFactor rScalingFactor
void NuTo::NodeBase::ScaleShapeFunctionMultiscalePeriodic(int rShapeFunction, double rScalingFactor)
{
	throw MechanicsException("[NuTo::NodeBase::ScaleShapeFunctionMultiscalePeriodic] Node of type " + GetNodeTypeStr() + " is not a multiscale node.");
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
				switch (this->GetNumVelocities())
				{
				case 1:
					this->GetVelocities1D(velocities);
					break;
				case 2:
					this->GetVelocities2D(velocities);
					break;
				case 3:
					this->GetVelocities3D(velocities);
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
				switch (this->GetNumAngularVelocities())
				{
				case 1:
					this->GetAngularVelocities2D(angularVelocities);
					break;
				case 3:
					this->GetAngularVelocities3D(angularVelocities);
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
				switch (this->GetNumAccelerations())
				{
				case 1:
					this->GetAccelerations1D(accelerations);
					break;
				case 2:
					this->GetAccelerations2D(accelerations);
					break;
				case 3:
					this->GetAccelerations3D(accelerations);
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
				switch (this->GetNumAngularAccelerations())
				{
				case 1:
					this->GetAngularAccelerations2D(angularAccelerations);
					break;
				case 3:
					this->GetAngularAccelerations3D(angularAccelerations);
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

