// $Id: $
#ifndef NodeCoordinatesTemperatures2D_H
#define NodeCoordinatesTemperatures2D_H

#include "nuto/mechanics/nodes/NodeCoordinates2D.h"
#include "nuto/mechanics/nodes/NodeTemperature.h"

namespace NuTo
{
//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... class for nodes having coordinates and Temperature
class NodeCoordinatesTemperature2D : public NodeCoordinates2D, public NodeTemperature
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:

    //! @brief constructor
    NodeCoordinatesTemperature2D() : NodeCoordinates2D (), NodeTemperature()
    {}

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeCoordinates2D)
           & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeTemperature);
    }
#endif  // ENABLE_SERIALIZATION
    //! @brief sets the global dofs
    //! @param rDOF current maximum DOF, this variable is increased within the routine
    virtual void SetGlobalDofs(int& rDOF)
    {
        NodeCoordinates2D::SetGlobalDofs(rDOF);
        NodeTemperature::SetGlobalDofs(rDOF);
    }

    //! @brief write dof values to the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void SetGlobalDofValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
    {
        NodeCoordinates2D::SetGlobalDofValues(rActiveDofValues, rDependentDofValues);
        NodeTemperature::SetGlobalDofValues(rActiveDofValues, rDependentDofValues);
    }

    //! @brief extract dof values from the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void GetGlobalDofValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
    {
        NodeCoordinates2D::GetGlobalDofValues(rActiveDofValues, rDependentDofValues);
        NodeTemperature::GetGlobalDofValues(rActiveDofValues, rDependentDofValues);
    }

    //! @brief write first time derivative of the dof values (e.g. velocities) to the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void SetGlobalDofFirstTimeDerivativeValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
    {
        NodeCoordinates2D::SetGlobalDofFirstTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
        NodeTemperature::SetGlobalDofFirstTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
    }

    //! @brief extract first time derivative of the dof values (e.g. velocities) from the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void GetGlobalDofFirstTimeDerivativeValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
    {
        NodeCoordinates2D::GetGlobalDofFirstTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
        NodeTemperature::GetGlobalDofFirstTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
    }

    //! @brief write second time derivative of the dof values (e.g. accelerations) to the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void SetGlobalDofSecondTimeDerivativeValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
    {
        NodeCoordinates2D::SetGlobalDofSecondTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
        NodeTemperature::SetGlobalDofSecondTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
    }

    //! @brief extract second time derivative of the dof values (e.g. accelerations) from the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void GetGlobalDofSecondTimeDerivativeValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
    {
        NodeCoordinates2D::GetGlobalDofSecondTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
        NodeTemperature::GetGlobalDofSecondTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
    }

    //! @brief renumber the global dofs according to predefined ordering
    //! @param rMappingInitialToNewOrdering ... mapping from initial ordering to the new ordering
    virtual void RenumberGlobalDofs(std::vector<int>& rMappingInitialToNewOrdering)
    {
        NodeCoordinates2D::RenumberGlobalDofs(rMappingInitialToNewOrdering);
        NodeTemperature::RenumberGlobalDofs(rMappingInitialToNewOrdering);
    }

    //! @brief returns the type of the node
    //! @return type
    virtual std::string GetNodeTypeStr()const
    {
    	return std::string("NodeCoordinatesTemperature2D");
    }
};
}

#endif //NodeCoordinatesTemperatures2D_H
