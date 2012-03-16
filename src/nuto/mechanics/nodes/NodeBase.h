// $Id$
#ifndef NODEBASE_H
#define NODEBASE_H

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif  // ENABLE_SERIALIZATION

#include <vector>
#include <boost/array.hpp>
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/nodes/NodeEnum.h"

#ifdef ENABLE_VISUALIZE
#include "nuto/visualize/VisualizeBase.h"
#include "nuto/visualize/VisualizeComponentBase.h"
#include "nuto/visualize/VisualizeUnstructuredGrid.h"
#include <boost/ptr_container/ptr_list.hpp>
#endif // ENABLE_VISUALIZE

namespace NuTo
{
template <class T> class FullMatrix;
class NodeDisplacementsMultiscale2D;
//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... standard abstract class for all nodes
class NodeBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
    //! @brief constructor
    NodeBase();

    //! @brief destructor
    virtual ~NodeBase(){};

    //! @brief assignment operator
    void operator=(NodeBase const& rOther)
    {
    }

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION

    //! @brief sets the global dofs
    //! @param rDOF current maximum DOF, this variable is increased within the routine
    virtual void SetGlobalDofs(int& rDOF)=0;

    //! @brief write dof values to the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void SetGlobalDofValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues) = 0;

    //! @brief extract dof values from the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void GetGlobalDofValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const = 0;

    //! @brief extract all dof numbers from the node (based on global dof number)
    virtual int* GetGlobalDofs();

    		//! @brief write first time derivative of the dof values (e.g. velocities) to the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void SetGlobalDofFirstTimeDerivativeValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues) = 0;

    //! @brief extract first time derivative of the dof values (e.g. velocities) from the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void GetGlobalDofFirstTimeDerivativeValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const = 0;

    //! @brief write second time derivative of the dof values (e.g. accelerations) to the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void SetGlobalDofSecondTimeDerivativeValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues) = 0;

    //! @brief extract second time derivative of the dof values (e.g. accelerations) from the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void GetGlobalDofSecondTimeDerivativeValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const = 0;

    //! @brief renumber the global dofs according to predefined ordering
    //! @param rMappingInitialToNewOrdering ... mapping from initial ordering to the new ordering
    virtual void RenumberGlobalDofs(std::vector<int>& rMappingInitialToNewOrdering) = 0;

    //! @brief returns the number of coordinates of the node
    //! @return number of coordinates
    virtual int GetNumCoordinates()const;

    //! @brief returns the coordinates of the node
    //! @return coordinates
    virtual void GetCoordinates1D(double rCoordinates[1])const;

    //! @brief set the coordinates
    //! @param rCoordinates  given coordinates
    virtual void SetCoordinates1D(const double rCoordinates[1]);

    //! @brief returns the coordinates of the node
    //! @return coordinates
    virtual void GetCoordinates2D(double rCoordinates[2])const;

    //! @brief set the coordinates
    //! @param rCoordinates  given coordinates
    virtual void SetCoordinates2D(const double rCoordinates[2]);

    //! @brief returns the coordinates of the node
    //! @return coordinates
    virtual void GetCoordinates3D(double rCoordinates[3])const;

    //! @brief set the coordinates
    //! @param rCoordinates  given coordinates
    virtual void SetCoordinates3D(const double rCoordinates[3]);

    //! @brief returns the number of coordinates of the node
    //! @return coordinates
    virtual double GetCoordinate(short rIndex)const;

    //! @brief returns the number of displacements of the node
    //! @return number of displacements
    virtual int GetNumDisplacements()const;

    //! @brief gives the global DOF of a displacement component
    //! @param rComponent component
    //! @return global DOF
    virtual int GetDofDisplacement(int rComponent)const;

    //! @brief returns the displacements of the node
    //! @return displacement
    virtual void GetDisplacements1D(double rCoordinates[1])const;

    //! @brief set the displacements
    //! @param rDisplacements  given displacements
    virtual void SetDisplacements1D(const double rDisplacements[1]);

    //! @brief returns the displacements of the node
    //! @return displacement
    virtual void GetDisplacements2D(double rCoordinates[2])const;

    //! @brief set the displacements
    //! @param rDisplacements  given displacements
    virtual void SetDisplacements2D(const double rDisplacements[2]);

    //! @brief returns the displacements of the node
    //! @return displacement
    virtual void GetDisplacements3D(double rCoordinates[3])const;

    //! @brief set the displacements
    //! @param rDisplacements  given displacements
    virtual void SetDisplacements3D(const double rDisplacements[3]);

    //! @brief returns the displacements of the node
    //! @return displacement
    virtual double GetDisplacement(short rIndex)const;

    //! @brief returns the number of displacements of the node
    //! @return number of displacements
    virtual int GetNumFineScaleDisplacements()const;

    //! @brief gives the global DOF of a displacement component
    //! @param rComponent component
    //! @return global DOF
    virtual int GetDofFineScaleDisplacement(int rComponent)const;

    //! @brief set the displacements
    //! @param rDisplacements  given displacements
    virtual void SetFineScaleDisplacements2D(const double rDisplacements[2]);

    //! @brief writes the displacements of a node to the prescribed pointer
    //! @param rDisplacements displacements
    virtual void GetFineScaleDisplacements2D(double rDisplacements[2])const;

    //! @brief returns the number of velocities of the node
    //! @return number of velocities
    virtual int GetNumVelocities()const;

    //! @brief gives the global DOF of a velocity component (Note: the velocity dof number is derived from the displacement dof)
    //! @param rComponent component
    //! @return global DOF
    virtual int GetDofVelocity(int rComponent)const;

    //! @brief returns the velocities of the node
    //! @param rVelocities ... velocities
    virtual void GetVelocities1D(double rVelocities[1])const;

    //! @brief set the velocities
    //! @param rVelocities  given velocities
    virtual void SetVelocities1D(const double rVelocities[1]);

    //! @brief returns the velocities of the node
    //! @param rVelocities ... velocities
    virtual void GetVelocities2D(double rVelocities[2])const;

    //! @brief set the velocities
    //! @param rVelocities  given velocities
    virtual void SetVelocities2D(const double rVelocities[2]);

    //! @brief returns the velocities of the node
    //! @param rVelocities ... velocities
    virtual void GetVelocities3D(double rVelocities[3])const;

    //! @brief set the velocities
    //! @param rVelocities  given velocities
    virtual void SetVelocities3D(const double rVelocities[3]);

    //! @brief returns the number of accelerations of the node
    //! @return number of accelerations
    virtual int GetNumAccelerations()const;

    //! @brief gives the global DOF of a acceleration component (Note: the acceleration dof number is derived from the displacement dof)
    //! @param rComponent component
    //! @return global DOF
    virtual int GetDofAcceleration(int rComponent)const;

    //! @brief returns the accelerations of the node
    //! @param rAccelerations ... accelerations
    virtual void GetAccelerations1D(double rAccelerations[1])const;

    //! @brief set the accelerations
    //! @param rVelocities  given accelerations
    virtual void SetAccelerations1D(const double rAccelerations[1]);

    //! @brief returns the accelerations of the node
    //! @param rAccelerations ... accelerations
    virtual void GetAccelerations2D(double rAccelerations[2])const;

    //! @brief set the accelerations
    //! @param rVelocities  given accelerations
    virtual void SetAccelerations2D(const double rAccelerations[2]);

    //! @brief returns the accelerations of the node
    //! @param rAccelerations ... accelerations
    virtual void GetAccelerations3D(double rAccelerations[3])const;

    //! @brief set the accelerations
    //! @param rVelocities  given accelerations
    virtual void SetAccelerations3D(const double rAccelerations[3]);

    //! @brief returns the number of Rotations of the node
    //! @return number of Rotations
    virtual int GetNumRotations()const;

    //! @brief gives the global DOF of a Rotation component
    //! @param rComponent component
    //! @return global DOF
    virtual int GetDofRotation(int rComponent)const;

    //! @brief returns the Rotations of the node
    //! @return Rotation
    virtual void GetRotations2D(double rCoordinates[2])const;

    //! @brief set the Rotations
    //! @param rRotations  given Rotations
    virtual void SetRotations2D(const double rRotations[2]);

    //! @brief returns the Rotations of the node
    //! @return Rotation
    virtual void GetRotations3D(double rCoordinates[3])const;

    //! @brief set the Rotations
    //! @param rRotations  given Rotations
    virtual void SetRotations3D(const double rRotations[3]);

    //! @brief returns the Rotations of the node
    //! @return Rotation
    virtual double GetRotation(short rIndex)const;

    //! @brief returns the number of radii
    //! @return number of radii
    virtual int GetNumRadius()const;

    //! @brief returns the radius of the node
    //! @param rRadius ... radius
    virtual void GetRadius(double rRadius[1])const;

    //! @brief set the radius
    //! @param rRadius  given radius
    virtual void SetRadius(const double rRadius[1]);

    //! @brief returns the number of temperatures of the node
    //! @return number of temperatures
    virtual int GetNumTemperatures()const;

    //! @brief gives the global DOF of a temperature component
    //! @param rComponent component
    //! @return global DOF
    virtual int GetDofTemperature(int rComponent)const;

    //! @brief returns the type of node as a string (all the data stored at the node)
    //! @return string
    virtual std::string GetNodeTypeStr()const=0;

    //! @brief returns the type of node as an enum (all the data stored at the node)
    //! @return enum
    virtual Node::eNodeType GetNodeType()const=0;

#ifdef ENABLE_VISUALIZE
    virtual void Visualize(VisualizeUnstructuredGrid& rVisualize, const boost::ptr_list<NuTo::VisualizeComponentBase>& rWhat) const;

#endif // ENABLE_VISUALIZE

    //! @brief set the shape functions based on the actual oscillations
    //! @parameter shape function number (0..2)
    virtual void SetShapeFunctionMultiscalePeriodic(int rShapeFunction);

    //! @brief returns the shape function for the periodic bc for the nodes
    virtual const boost::array<double,3>& GetShapeFunctionMultiscalePeriodicX()const;

    //! @brief returns the shape function for the periodic bc for the nodes
    virtual const boost::array<double,3>& GetShapeFunctionMultiscalePeriodicY() const;

    //! @brief scales the shape functions
    //! @parameter rShapeFunction  (1..3 corresponding to macro strains exx, eyy, and gxy)
    //! @parameter rScalingFactor rScalingFactor
    virtual void ScaleShapeFunctionMultiscalePeriodic(int rShapeFunction, double rScalingFactor);

    //! @brief only relevant for multiscale nodes, where they are either in the homogeneous (false) or cracked domain (true)
    virtual bool IsInCrackedDomain()const
    {
    	throw MechanicsException("[NuTo::NodeBase::IsInCrackedDomain] not implemented for this type of nodes.");
    }

    virtual NodeDisplacementsMultiscale2D* AsNodeDisplacementsMultiscale2D()
    {
    	throw MechanicsException("[NuTo::NodeBase::AsNodeDisplacementsMultiscale2D] conversion can't be performed, types do not match.");
    }

    virtual const NodeDisplacementsMultiscale2D* AsNodeDisplacementsMultiscale2D()const
    {
    	throw MechanicsException("[NuTo::NodeBase::AsNodeDisplacementsMultiscale2D] conversion can't be performed, types do not match.");
    }

    //! @brief clones (copies) the node with all its data, it's supposed to be a new node, so be careful with ptr
    virtual NodeBase* Clone()const=0;

protected:
    //the base class of the nodes must not contain any data

};
}//namespace NuTo
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::NodeBase)
#endif // ENABLE_SERIALIZATION
#endif //NODEBASE_H

