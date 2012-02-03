// $Id:$
#ifndef NodeCoordinatesDisplacementsRotationsRadius_2d_H
#define NodeCoordinatesDisplacementsRotationsRadius_2d_H

#include "nuto/mechanics/nodes/NodeCoordinates2D.h"
#include "nuto/mechanics/nodes/NodeDisplacements2D.h"
#include "nuto/mechanics/nodes/NodeRotations2D.h"
#include "nuto/mechanics/nodes/NodeRadius.h"

namespace NuTo
{
//! @author Jörg F. Unger, NU
//! @date January 2012
//! @brief ... class for nodes having coordinates, displacements and a radius
class NodeCoordinatesDisplacementsRotationsRadius2D : public NodeCoordinates2D, public NodeDisplacements2D, public NodeRotations2D, public NodeRadius
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:

    //! @brief constructor
    NodeCoordinatesDisplacementsRotationsRadius2D();

    //! @brief copy constructor
    NodeCoordinatesDisplacementsRotationsRadius2D(NodeCoordinatesDisplacementsRotationsRadius2D const& rOther)
    {
        (*this) = rOther;
    }

    //! @brief clones (copies) the node with all its data, it's supposed to be a new node, so be careful with ptr
    NodeCoordinatesDisplacementsRotationsRadius2D* Clone()const
    {
    	return new NodeCoordinatesDisplacementsRotationsRadius2D(*this);
    }

    //! @brief assignment operator
    NodeCoordinatesDisplacementsRotationsRadius2D& operator= (NodeCoordinatesDisplacementsRotationsRadius2D const& rOther);

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
        NodeCoordinates2D::SetGlobalDofs(rDOF);
        NodeDisplacements2D::SetGlobalDofs(rDOF);
        NodeRotations2D::SetGlobalDofs(rDOF);
        NodeRadius::SetGlobalDofs(rDOF);
    }

    //! @brief write dof values to the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void SetGlobalDofValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
    {
        NodeCoordinates2D::SetGlobalDofValues(rActiveDofValues, rDependentDofValues);
        NodeDisplacements2D::SetGlobalDofValues(rActiveDofValues, rDependentDofValues);
        NodeRotations2D::SetGlobalDofValues(rActiveDofValues, rDependentDofValues);
        NodeRadius::SetGlobalDofValues(rActiveDofValues, rDependentDofValues);
    }

    //! @brief extract dof values from the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void GetGlobalDofValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
    {
        NodeCoordinates2D::GetGlobalDofValues(rActiveDofValues, rDependentDofValues);
        NodeDisplacements2D::GetGlobalDofValues(rActiveDofValues, rDependentDofValues);
        NodeRotations2D::GetGlobalDofValues(rActiveDofValues, rDependentDofValues);
        NodeRadius::GetGlobalDofValues(rActiveDofValues, rDependentDofValues);
    }

    //! @brief write first time derivative of the dof values (e.g. velocities) to the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void SetGlobalDofFirstTimeDerivativeValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
    {
        NodeCoordinates2D::SetGlobalDofFirstTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
        NodeDisplacements2D::SetGlobalDofFirstTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
        NodeRotations2D::SetGlobalDofFirstTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
        NodeRadius::SetGlobalDofFirstTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
    }

    //! @brief extract first time derivative of the dof values (e.g. velocities) from the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void GetGlobalDofFirstTimeDerivativeValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
    {
        NodeCoordinates2D::GetGlobalDofFirstTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
        NodeDisplacements2D::GetGlobalDofFirstTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
        NodeRotations2D::GetGlobalDofFirstTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
        NodeRadius::GetGlobalDofFirstTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
    }

    //! @brief write second time derivative of the dof values (e.g. accelerations) to the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void SetGlobalDofSecondTimeDerivativeValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
    {
        NodeCoordinates2D::SetGlobalDofSecondTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
        NodeDisplacements2D::SetGlobalDofSecondTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
        NodeRotations2D::SetGlobalDofSecondTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
        NodeRadius::SetGlobalDofSecondTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
    }

    //! @brief extract second time derivative of the dof values (e.g. accelerations) from the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void GetGlobalDofSecondTimeDerivativeValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
    {
        NodeCoordinates2D::GetGlobalDofSecondTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
        NodeDisplacements2D::GetGlobalDofSecondTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
        NodeRotations2D::GetGlobalDofSecondTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
        NodeRadius::GetGlobalDofSecondTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
    }

    //! @brief renumber the global dofs according to predefined ordering
    //! @param rMappingInitialToNewOrdering ... mapping from initial ordering to the new ordering
    virtual void RenumberGlobalDofs(std::vector<int>& rMappingInitialToNewOrdering)
    {
        NodeCoordinates2D::RenumberGlobalDofs(rMappingInitialToNewOrdering);
        NodeDisplacements2D::RenumberGlobalDofs(rMappingInitialToNewOrdering);
        NodeRotations2D::RenumberGlobalDofs(rMappingInitialToNewOrdering);
        NodeRadius::RenumberGlobalDofs(rMappingInitialToNewOrdering);
    }

    //! @brief returns the type of the node
    //! @return type
    virtual std::string GetNodeTypeStr()const
    {
    	return std::string("NodeCoordinatesDisplacementsRotationsRadius2D");
    }

    //! @brief returns the type of node as an enum (all the data stored at the node)
    //! @return enum
    virtual Node::eNodeType GetNodeType()const
    {
        return Node::NodeCoordinatesDisplacementsRotationsRadius2D;
    }
};
}
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::NodeCoordinatesDisplacementsRotationsRadius2D)
namespace boost{
template<>
struct is_virtual_base_of<NuTo::NodeCoordinates2D, NuTo::NodeCoordinatesDisplacementsRotationsRadius2D>: public mpl::true_ {};
template<>
struct is_virtual_base_of<NuTo::NodeDisplacements2D, NuTo::NodeCoordinatesDisplacementsRotationsRadius2D>: public mpl::true_ {};
template<>
struct is_virtual_base_of<NuTo::NodeRotations2D, NuTo::NodeCoordinatesDisplacementsRotationsRadius2D>: public mpl::true_ {};
template<>
struct is_virtual_base_of<NuTo::NodeRadius, NuTo::NodeCoordinatesDisplacementsRotationsRadius2D>: public mpl::true_ {};
}
#endif // ENABLE_SERIALIZATION

#endif //NodeCoordinatesDisplacementsRotationsRadius_2d_H
