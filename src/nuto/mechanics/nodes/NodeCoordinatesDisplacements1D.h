// $Id$
#ifndef NodeCoordinatesDisplacements_1d_H
#define NodeCoordinatesDisplacements_1d_H

#include "nuto/mechanics/nodes/NodeCoordinates1D.h"
#include "nuto/mechanics/nodes/NodeDisplacements1D.h"

namespace NuTo
{
//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... class for nodes having coordinates and displacements
class NodeCoordinatesDisplacements1D : public NodeCoordinates1D, public NodeDisplacements1D
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:

    //! @brief constructor
    NodeCoordinatesDisplacements1D();

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
        NodeCoordinates1D::SetGlobalDofs(rDOF);
        NodeDisplacements1D::SetGlobalDofs(rDOF);
    }

    //! @brief write dof values to the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void SetGlobalDofValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
    {
        NodeCoordinates1D::SetGlobalDofValues(rActiveDofValues, rDependentDofValues);
        NodeDisplacements1D::SetGlobalDofValues(rActiveDofValues, rDependentDofValues);
    }

    //! @brief extract dof values from the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void GetGlobalDofValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
    {
        NodeCoordinates1D::GetGlobalDofValues(rActiveDofValues, rDependentDofValues);
        NodeDisplacements1D::GetGlobalDofValues(rActiveDofValues, rDependentDofValues);
    }

    //! @brief write first time derivative of the dof values (e.g. velocities) to the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void SetGlobalDofFirstTimeDerivativeValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
    {
        NodeCoordinates1D::SetGlobalDofFirstTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
        NodeDisplacements1D::SetGlobalDofFirstTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
    }

    //! @brief extract first time derivative of the dof values (e.g. velocities) from the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void GetGlobalDofFirstTimeDerivativeValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
    {
        NodeCoordinates1D::GetGlobalDofFirstTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
        NodeDisplacements1D::GetGlobalDofFirstTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
    }

    //! @brief write second time derivative of the dof values (e.g. accelerations) to the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void SetGlobalDofSecondTimeDerivativeValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
    {
        NodeCoordinates1D::SetGlobalDofSecondTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
        NodeDisplacements1D::SetGlobalDofSecondTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
    }

    //! @brief extract second time derivative of the dof values (e.g. accelerations) from the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void GetGlobalDofSecondTimeDerivativeValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
    {
        NodeCoordinates1D::GetGlobalDofSecondTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
        NodeDisplacements1D::GetGlobalDofSecondTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
    }

    //! @brief renumber the global dofs according to predefined ordering
    //! @param rMappingInitialToNewOrdering ... mapping from initial ordering to the new ordering
    virtual void RenumberGlobalDofs(std::vector<int>& rMappingInitialToNewOrdering)
    {
        NodeCoordinates1D::RenumberGlobalDofs(rMappingInitialToNewOrdering);
        NodeDisplacements1D::RenumberGlobalDofs(rMappingInitialToNewOrdering);
    }

    //! @brief returns the type of the node
    //! @return type
    virtual std::string GetNodeTypeStr()const
    {
    	return std::string("NodeCoordinatesDisplacements1D");
    }
};
}
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::NodeCoordinatesDisplacements1D)
namespace boost{
template<>
struct is_virtual_base_of<NuTo::NodeCoordinates1D, NuTo::NodeCoordinatesDisplacements1D>: public mpl::true_ {};
template<>
struct is_virtual_base_of<NuTo::NodeDisplacements1D, NuTo::NodeCoordinatesDisplacements1D>: public mpl::true_ {};
}
#endif // ENABLE_SERIALIZATION

#endif //NodeCoordinatesDisplacements_H
