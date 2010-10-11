// $Id$
#ifndef NodeCoordinatesDisplacementsRotations3D_H
#define NodeCoordinatesDisplacementsRotations3D_H

#include "nuto/mechanics/nodes/NodeCoordinates3D.h"
#include "nuto/mechanics/nodes/NodeDisplacements3D.h"
#include "nuto/mechanics/nodes/NodeRotations3D.h"

namespace NuTo
{
//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... class for nodes having coordinates, displacements and rotations
class NodeCoordinatesDisplacementsRotations3D : public NodeCoordinates3D,public NodeDisplacements3D, public NodeRotations3D
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:

    //! @brief constructor
    NodeCoordinatesDisplacementsRotations3D() : NodeCoordinates3D (),
            NodeDisplacements3D(),
            NodeRotations3D()
    {}

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeCoordinates3D)
           & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeDisplacements3D)
           & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeRotations3D);
    }
#endif  // ENABLE_SERIALIZATION
    //! @brief sets the global dofs
    //! @param rDOF current maximum DOF, this variable is increased within the routine
    virtual void SetGlobalDofs(int& rDOF)
    {
        NodeCoordinates3D::SetGlobalDofs(rDOF);
        NodeDisplacements3D::SetGlobalDofs(rDOF);
        NodeRotations3D::SetGlobalDofs(rDOF);
    }

    //! @brief write dof values to the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void SetGlobalDofValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
    {
        NodeCoordinates3D::SetGlobalDofValues(rActiveDofValues, rDependentDofValues);
        NodeDisplacements3D::SetGlobalDofValues(rActiveDofValues, rDependentDofValues);
        NodeRotations3D::SetGlobalDofValues(rActiveDofValues, rDependentDofValues);
    }

    //! @brief extract dof values from the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void GetGlobalDofValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
    {
        NodeCoordinates3D::GetGlobalDofValues(rActiveDofValues, rDependentDofValues);
        NodeDisplacements3D::GetGlobalDofValues(rActiveDofValues, rDependentDofValues);
        NodeRotations3D::GetGlobalDofValues(rActiveDofValues, rDependentDofValues);
    }

    //! @brief write first time derivative of the dof values (e.g. velocities) to the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void SetGlobalDofFirstTimeDerivativeValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
    {
        NodeCoordinates3D::SetGlobalDofFirstTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
        NodeDisplacements3D::SetGlobalDofFirstTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
        NodeRotations3D::SetGlobalDofFirstTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
    }

    //! @brief extract first time derivative of the dof values (e.g. velocities) from the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void GetGlobalDofFirstTimeDerivativeValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
    {
        NodeCoordinates3D::GetGlobalDofFirstTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
        NodeDisplacements3D::GetGlobalDofFirstTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
        NodeRotations3D::GetGlobalDofFirstTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
    }

    //! @brief write second time derivative of the dof values (e.g. accelerations) to the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void SetGlobalDofSecondTimeDerivativeValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
    {
        NodeCoordinates3D::SetGlobalDofSecondTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
        NodeDisplacements3D::SetGlobalDofSecondTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
        NodeRotations3D::SetGlobalDofSecondTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
    }

    //! @brief extract second time derivative of the dof values (e.g. accelerations) from the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void GetGlobalDofSecondTimeDerivativeValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
    {
        NodeCoordinates3D::GetGlobalDofSecondTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
        NodeDisplacements3D::GetGlobalDofSecondTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
        NodeRotations3D::GetGlobalDofSecondTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
    }

    //! @brief renumber the global dofs according to predefined ordering
    //! @param rMappingInitialToNewOrdering ... mapping from initial ordering to the new ordering
    virtual void RenumberGlobalDofs(std::vector<int>& rMappingInitialToNewOrdering)
    {
        NodeCoordinates3D::RenumberGlobalDofs(rMappingInitialToNewOrdering);
        NodeDisplacements3D::RenumberGlobalDofs(rMappingInitialToNewOrdering);
        NodeRotations3D::RenumberGlobalDofs(rMappingInitialToNewOrdering);
    }

    //! @brief returns the type of the node
    //! @return type
    virtual std::string GetNodeTypeStr()const
    {
    	return std::string("NodeCoordinatesDisplacementsRotations3D");
    }
};
}

#endif //NodeCoordinatesDisplacementsRotations3D_H
