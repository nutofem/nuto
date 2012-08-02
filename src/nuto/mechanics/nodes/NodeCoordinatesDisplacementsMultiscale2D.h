// $Id: NodeCoordinatesDisplacementsMultiscale2D.h 331 2010-10-06 09:32:11Z arnold2 $
#ifndef NodeCoordinatesDisplacementsMultiscale2D_H
#define NodeCoordinatesDisplacementsMultiscale2Dd_H

#include "nuto/mechanics/nodes/NodeCoordinates_Def.h"
#include "nuto/mechanics/nodes/NodeDisplacementsMultiscale2D.h"

namespace NuTo
{
class NodeCoordinatesDisplacements2D;
//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... class for nodes having coordinates and displacements
class NodeCoordinatesDisplacementsMultiscale2D : public NodeCoordinates<2>, public NodeDisplacementsMultiscale2D
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:

    //! @brief constructor
    NodeCoordinatesDisplacementsMultiscale2D(NuTo::StructureMultiscale* rStructureMultiscale, bool rCrackedDomain);

    //! @brief copy constructor
    NodeCoordinatesDisplacementsMultiscale2D(NodeCoordinatesDisplacementsMultiscale2D const& rOther)
    {
        (*this) = rOther;
    }

    //! @brief clones (copies) the node with all its data, it's supposed to be a new node, so be careful with ptr
    NodeCoordinatesDisplacementsMultiscale2D* Clone()const
    {
    	return new NodeCoordinatesDisplacementsMultiscale2D(*this);
    }

    //! @brief assignment operator
    NodeCoordinatesDisplacementsMultiscale2D& operator= (NodeCoordinatesDisplacementsMultiscale2D const& rOther);

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
        NodeDisplacementsMultiscale2D::SetGlobalDofs(rDOF);
    }

    //! @brief write dof values to the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void SetGlobalDofValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
    {
        NodeDisplacementsMultiscale2D::SetGlobalDofValues(rActiveDofValues, rDependentDofValues);
    }

    //! @brief extract dof values from the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void GetGlobalDofValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
    {
        NodeDisplacementsMultiscale2D::GetGlobalDofValues(rActiveDofValues, rDependentDofValues);
    }

    //! @brief write first time derivative of the dof values (e.g. velocities) to the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void SetGlobalDofFirstTimeDerivativeValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
    {
        NodeDisplacementsMultiscale2D::SetGlobalDofFirstTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
    }

    //! @brief extract first time derivative of the dof values (e.g. velocities) from the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void GetGlobalDofFirstTimeDerivativeValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
    {
        NodeDisplacementsMultiscale2D::GetGlobalDofFirstTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
    }

    //! @brief write second time derivative of the dof values (e.g. accelerations) to the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void SetGlobalDofSecondTimeDerivativeValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
    {
        NodeDisplacementsMultiscale2D::SetGlobalDofSecondTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
    }

    //! @brief extract second time derivative of the dof values (e.g. accelerations) from the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void GetGlobalDofSecondTimeDerivativeValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
    {
        NodeDisplacementsMultiscale2D::GetGlobalDofSecondTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
    }

    //! @brief renumber the global dofs according to predefined ordering
    //! @param rMappingInitialToNewOrdering ... mapping from initial ordering to the new ordering
    virtual void RenumberGlobalDofs(std::vector<int>& rMappingInitialToNewOrdering)
    {
        NodeDisplacementsMultiscale2D::RenumberGlobalDofs(rMappingInitialToNewOrdering);
    }

    //! @brief returns the type of the node
    //! @return type
    virtual std::string GetNodeTypeStr()const
    {
    	return std::string("NodeCoordinatesDisplacementsMultiscale2D");
    }

/*    //! @brief returns the type of node as an enum (all the data stored at the node)
    //! @return enum
    virtual Node::eNodeType GetNodeType()const
    {
        return Node::NodeCoordinatesDisplacementsMultiscale2D;
    }
*/
protected:
    //! @brief constructor
    NodeCoordinatesDisplacementsMultiscale2D(){};
};
}
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::NodeCoordinatesDisplacementsMultiscale2D)
namespace boost{
template<>
struct is_virtual_base_of<NuTo::NodeCoordinates<2>, NuTo::NodeCoordinatesDisplacementsMultiscale2D>: public mpl::true_ {};
template<>
struct is_virtual_base_of<NuTo::NodeDisplacementsMultiscale2D, NuTo::NodeCoordinatesDisplacementsMultiscale2D>: public mpl::true_ {};
}
#endif // ENABLE_SERIALIZATION

#endif //NodeCoordinatesDisplacementsMultiscale2D_H
