// $Id$

#ifndef NODE_COORDINATES_DISPLACEMENTS_VELOCITIES_ACCELERATIONS_ROTATIOND_ANGULAR_VELOCITIES_ANGULAR_ACCELERATIONS_RADIUS2D_H_
#define NODE_COORDINATES_DISPLACEMENTS_VELOCITIES_ACCELERATIONS_ROTATIOND_ANGULAR_VELOCITIES_ANGULAR_ACCELERATIONS_RADIUS2D_H_

#include "nuto/mechanics/nodes/NodeCoordinates2D.h"
#include "nuto/mechanics/nodes/NodeDisplacements2D.h"
#include "nuto/mechanics/nodes/NodeVelocities2D.h"
#include "nuto/mechanics/nodes/NodeAccelerations2D.h"
#include "nuto/mechanics/nodes/NodeRotations2D.h"
#include "nuto/mechanics/nodes/NodeAngularVelocities2D.h"
#include "nuto/mechanics/nodes/NodeAngularAccelerations2D.h"
#include "nuto/mechanics/nodes/NodeRadius.h"

namespace NuTo
{
//! @author Joerg F. Unger
//! @date March 2012
//! @brief ... class for nodes having coordinates, displacements, velocities and accelerations
class NodeCoordinatesDisplacementsVelocitiesAccelerationsRotationsAngularVelocitiesAngularAccelerationsRadius2D : public virtual NuTo::NodeCoordinates2D,
    public NuTo::NodeDisplacements2D, public NuTo::NodeVelocities2D, public NuTo::NodeAccelerations2D,
    public NuTo::NodeRotations2D, public NuTo::NodeAngularVelocities2D, public NuTo::NodeAngularAccelerations2D, public NuTo::NodeRadius
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
    //! @brief constructor
    NodeCoordinatesDisplacementsVelocitiesAccelerationsRotationsAngularVelocitiesAngularAccelerationsRadius2D();

    //! @brief sets the global dofs
    //! @param rDOF current maximum DOF, this variable is increased within the routine
    virtual void SetGlobalDofs(int& rDOF);

    //! @brief write dof values to the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void SetGlobalDofValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues);

    //! @brief extract dof values from the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void GetGlobalDofValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const;

    //! @brief write first time derivative of the dof values (e.g. velocities) to the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void SetGlobalDofFirstTimeDerivativeValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues);

    //! @brief extract first time derivative of the dof values (e.g. velocities) from the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void GetGlobalDofFirstTimeDerivativeValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const;

    //! @brief write second time derivative of the dof values (e.g. accelerations) to the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void SetGlobalDofSecondTimeDerivativeValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues);

    //! @brief extract second time derivative of the dof values (e.g. accelerations) from the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void GetGlobalDofSecondTimeDerivativeValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const;

    //! @brief renumber the global dofs according to predefined ordering
    //! @param rMappingInitialToNewOrdering ... mapping from initial ordering to the new ordering
    virtual void RenumberGlobalDofs(std::vector<int>& rMappingInitialToNewOrdering);

    //! @brief returns the type of the node
    //! @return type
    virtual std::string GetNodeTypeStr()const;

    //! @brief returns the type of node as an enum (all the data stored at the node)
    //! @return enum
    virtual Node::eNodeType GetNodeType()const;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeCoordinates2D)
           & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeDisplacements2D)
           & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeVelocities2D)
           & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeAccelerations2D)
           & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeRotations2D)
           & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeAngularVelocities2D)
           & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeAngularAccelerations2D);
           & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeRadius);
    }
#endif  // ENABLE_SERIALIZATION
};

}
#endif // NODE_COORDINATES_DISPLACEMENTS_VELOCITIES_ACCELERATIONS_ROTATIOND_ANGULAR_VELOCITIES_ANGULAR_ACCELERATIONS_RADIUS2D_H_
