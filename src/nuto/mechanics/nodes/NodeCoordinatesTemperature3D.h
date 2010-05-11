// $Id: $
#ifndef NodeCoordinatesTemperatures3D_H
#define NodeCoordinatesTemperatures3D_H

#include "nuto/mechanics/nodes/NodeCoordinates3D.h"
#include "nuto/mechanics/nodes/NodeTemperature.h"

namespace NuTo
{
//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... class for nodes having coordinates and Temperature
class NodeCoordinatesTemperature3D : public  NodeCoordinates3D, public NodeTemperature
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:

    //! @brief constructor
    NodeCoordinatesTemperature3D() : NodeCoordinates3D (), NodeTemperature()
    {}

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeCoordinates3D)
           & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeTemperature);
    }
#endif  // ENABLE_SERIALIZATION
    //! @brief sets the global dofs
    //! @param rDOF current maximum DOF, this variable is increased within the routine
    virtual void SetGlobalDofs(int& rDOF)
    {
        NodeCoordinates3D::SetGlobalDofs(rDOF);
        NodeTemperature::SetGlobalDofs(rDOF);
    }

    //! @brief write dof values to the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void SetGlobalDofValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
    {
        NodeCoordinates3D::SetGlobalDofValues(rActiveDofValues, rDependentDofValues);
        NodeTemperature::SetGlobalDofValues(rActiveDofValues, rDependentDofValues);
    }

    //! @brief extract dof values from the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void GetGlobalDofValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
    {
        NodeCoordinates3D::GetGlobalDofValues(rActiveDofValues, rDependentDofValues);
        NodeTemperature::GetGlobalDofValues(rActiveDofValues, rDependentDofValues);
    }

    //! @brief renumber the global dofs according to predefined ordering
    //! @param rMappingInitialToNewOrdering ... mapping from initial ordering to the new ordering
    virtual void RenumberGlobalDofs(std::vector<int>& rMappingInitialToNewOrdering)
    {
        NodeCoordinates3D::RenumberGlobalDofs(rMappingInitialToNewOrdering);
        NodeTemperature::RenumberGlobalDofs(rMappingInitialToNewOrdering);
    }

    //! @brief returns the type of the node
    //! @return type
    virtual std::string GetNodeTypeStr()const
    {
    	return std::string("NodeCoordinatesTemperature3D");
    }
};
}

#endif //NodeCoordinatesTemperatures3D_H
