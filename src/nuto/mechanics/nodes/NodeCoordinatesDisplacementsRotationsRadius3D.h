// $Id:$
#ifndef NodeCoordinatesDisplacementsRotationsRadius_3d_H
#define NodeCoordinatesDisplacementsRotationsRadius_3d_H

#include "nuto/mechanics/nodes/NodeCoordinates3D.h"
#include "nuto/mechanics/nodes/NodeDisplacements3D.h"
#include "nuto/mechanics/nodes/NodeRotations3D.h"
#include "nuto/mechanics/nodes/NodeRadius.h"

namespace NuTo
{
//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... class for nodes having coordinates, displacements and a radius
class NodeCoordinatesDisplacementsRotationsRadius3D : public  NodeCoordinates3D,
    public NodeDisplacements3D,
    public NodeRotations3D,
    public NodeRadius
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
    //! @brief constructor
    NodeCoordinatesDisplacementsRotationsRadius3D();

    //! @brief copy constructor
    NodeCoordinatesDisplacementsRotationsRadius3D(NodeCoordinatesDisplacementsRotationsRadius3D const& rOther)
    {
        (*this) = rOther;
    }

    //! @brief clones (copies) the node with all its data, it's supposed to be a new node, so be careful with ptr
    NodeCoordinatesDisplacementsRotationsRadius3D* Clone()const
    {
    	return new NodeCoordinatesDisplacementsRotationsRadius3D(*this);
    }

    //! @brief assignment operator
    NodeCoordinatesDisplacementsRotationsRadius3D& operator= (NodeCoordinatesDisplacementsRotationsRadius3D const& rOther);


#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif  // ENABLE_SERIALIZATION

    //! @brief sets the global dofs
    //! @param rDOF current maximum DOF, this variable is increased within the routine
    virtual void SetGlobalDofs(int& rDOF)
    {
        NodeCoordinates3D::SetGlobalDofs(rDOF);
        NodeDisplacements3D::SetGlobalDofs(rDOF);
        NodeRotations3D::SetGlobalDofs(rDOF);
        NodeRadius::SetGlobalDofs(rDOF);
    }

    //! @brief write dof values to the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void SetGlobalDofValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
    {
        NodeCoordinates3D::SetGlobalDofValues(rActiveDofValues, rDependentDofValues);
        NodeDisplacements3D::SetGlobalDofValues(rActiveDofValues, rDependentDofValues);
        NodeRotations3D::SetGlobalDofValues(rActiveDofValues, rDependentDofValues);
        NodeRadius::SetGlobalDofValues(rActiveDofValues, rDependentDofValues);
    }

    //! @brief extract dof values from the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void GetGlobalDofValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
    {
        NodeCoordinates3D::GetGlobalDofValues(rActiveDofValues, rDependentDofValues);
        NodeDisplacements3D::GetGlobalDofValues(rActiveDofValues, rDependentDofValues);
        NodeRotations3D::GetGlobalDofValues(rActiveDofValues, rDependentDofValues);
        NodeRadius::GetGlobalDofValues(rActiveDofValues, rDependentDofValues);
    }

    //! @brief write first time derivative of the dof values (e.g. velocities) to the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void SetGlobalDofFirstTimeDerivativeValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
    {
        NodeCoordinates3D::SetGlobalDofFirstTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
        NodeDisplacements3D::SetGlobalDofFirstTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
        NodeRotations3D::SetGlobalDofFirstTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
        NodeRadius::SetGlobalDofFirstTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
    }

    //! @brief extract first time derivative of the dof values (e.g. velocities) from the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void GetGlobalDofFirstTimeDerivativeValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
    {
        NodeCoordinates3D::GetGlobalDofFirstTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
        NodeDisplacements3D::GetGlobalDofFirstTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
        NodeRotations3D::GetGlobalDofFirstTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
        NodeRadius::GetGlobalDofFirstTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
    }

    //! @brief write second time derivative of the dof values (e.g. accelerations) to the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void SetGlobalDofSecondTimeDerivativeValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
    {
        NodeCoordinates3D::SetGlobalDofSecondTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
        NodeDisplacements3D::SetGlobalDofSecondTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
        NodeRotations3D::SetGlobalDofSecondTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
        NodeRadius::SetGlobalDofSecondTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
    }

    //! @brief extract second time derivative of the dof values (e.g. accelerations) from the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void GetGlobalDofSecondTimeDerivativeValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
    {
        NodeCoordinates3D::GetGlobalDofSecondTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
        NodeDisplacements3D::GetGlobalDofSecondTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
        NodeRotations3D::GetGlobalDofSecondTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
        NodeRadius::GetGlobalDofSecondTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
    }

    //! @brief renumber the global dofs according to predefined ordering
    //! @param rMappingInitialToNewOrdering ... mapping from initial ordering to the new ordering
    virtual void RenumberGlobalDofs(std::vector<int>& rMappingInitialToNewOrdering)
    {
        NodeCoordinates3D::RenumberGlobalDofs(rMappingInitialToNewOrdering);
        NodeDisplacements3D::RenumberGlobalDofs(rMappingInitialToNewOrdering);
        NodeRotations3D::RenumberGlobalDofs(rMappingInitialToNewOrdering);
        NodeRadius::RenumberGlobalDofs(rMappingInitialToNewOrdering);
    }

    //! @brief returns the type of the node
    //! @return type
    virtual std::string GetNodeTypeStr()const
    {
    	return std::string("NodeCoordinatesDisplacementsRotationsRadius3D");
    }

    //! @brief returns the type of node as an enum (all the data stored at the node)
    //! @return enum
    virtual Node::eNodeType GetNodeType()const
    {
        return Node::NodeCoordinatesDisplacementsRotationsRadius3D;
    }
};
}
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::NodeCoordinatesDisplacementsRotationsRadius3D)
namespace boost{
template<>
struct is_virtual_base_of<NuTo::NodeCoordinates3D, NuTo::NodeCoordinatesDisplacementsRotationsRadius3D>: public mpl::true_ {};
template<>
struct is_virtual_base_of<NuTo::NodeDisplacements3D, NuTo::NodeCoordinatesDisplacementsRotationsRadius3D>: public mpl::true_ {};
template<>
struct is_virtual_base_of<NuTo::NodeRotations3D, NuTo::NodeCoordinatesDisplacementsRotationsRadius3D>: public mpl::true_ {};
template<>
struct is_virtual_base_of<NuTo::NodeRadius, NuTo::NodeCoordinatesDisplacementsRotationsRadius3D>: public mpl::true_ {};
}
#endif// ENABLE_SERIALIZATION

#endif //NodeCoordinatesDisplacementsRotationsRadius_3d_H
