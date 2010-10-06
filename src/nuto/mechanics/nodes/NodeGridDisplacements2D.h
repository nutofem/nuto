// $Id$
#ifndef NodeGridDisplacements_2d_H
#define NodeGridDisplacements_2d_H

#include "nuto/mechanics/nodes/NodeGrid2D.h"
#include "nuto/mechanics/nodes/NodeDisplacements2D.h"

namespace NuTo
{
//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... class for grid nodes with displacements
class NodeGridDisplacements2D : public NodeGrid2D, public NodeDisplacements2D
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:

    //! @brief constructor
    NodeGridDisplacements2D(int rNodeId) : NodeGrid2D (rNodeId), NodeDisplacements2D()
    {}

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeGrid2D)
           & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeDisplacements2D);
    }
#endif  // ENABLE_SERIALIZATION

    //! @brief sets the global dofs
    //! @param rDOF current maximum DOF, this variable is increased within the routine
    virtual void SetGlobalDofs(int& rDOF)
    {
        NodeGrid2D::SetGlobalDofs(rDOF);
        NodeDisplacements2D::SetGlobalDofs(rDOF);
    }

    //! @brief write dof values to the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void SetGlobalDofValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
    {
        NodeGrid2D::SetGlobalDofValues(rActiveDofValues, rDependentDofValues);
        NodeDisplacements2D::SetGlobalDofValues(rActiveDofValues, rDependentDofValues);
    }

    //! @brief extract dof values from the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void GetGlobalDofValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
    {
        NodeGrid2D::GetGlobalDofValues(rActiveDofValues, rDependentDofValues);
        NodeDisplacements2D::GetGlobalDofValues(rActiveDofValues, rDependentDofValues);
    }

    //! @brief write first time derivative of the dof values (e.g. velocities) to the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void SetGlobalDofFirstTimeDerivativeValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
    {
    	NodeGrid2D::SetGlobalDofFirstTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
        NodeDisplacements2D::SetGlobalDofFirstTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
    }

    //! @brief extract first time derivative of the dof values (e.g. velocities) from the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void GetGlobalDofFirstTimeDerivativeValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
    {
    	NodeGrid2D::GetGlobalDofFirstTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
        NodeDisplacements2D::GetGlobalDofFirstTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
    }

    //! @brief write second time derivative of the dof values (e.g. accelerations) to the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void SetGlobalDofSecondTimeDerivativeValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
    {
    	NodeGrid2D::SetGlobalDofSecondTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
        NodeDisplacements2D::SetGlobalDofSecondTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
    }

    //! @brief extract second time derivative of the dof values (e.g. accelerations) from the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void GetGlobalDofSecondTimeDerivativeValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
    {
    	NodeGrid2D::GetGlobalDofSecondTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
        NodeDisplacements2D::GetGlobalDofSecondTimeDerivativeValues(rActiveDofValues, rDependentDofValues);
    }

    //! @brief renumber the global dofs according to predefined ordering
    //! @param rMappingInitialToNewOrdering ... mapping from initial ordering to the new ordering
    virtual void RenumberGlobalDofs(std::vector<int>& rMappingInitialToNewOrdering)
    {
        NodeGrid2D::RenumberGlobalDofs(rMappingInitialToNewOrdering);
        NodeDisplacements2D::RenumberGlobalDofs(rMappingInitialToNewOrdering);
    }

    //! @brief returns the type of the node
    //! @return type
    virtual std::string GetNodeTypeStr()const
    {
    	return std::string("NodeGridDisplacements2D");
    }
};
}

#endif //NodeGridDisplacements_2d_H
