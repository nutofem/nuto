// $Id: $
#ifndef NodeGridDisplacements_1d_H
#define NodeGridDisplacements_1d_H

#include "nuto/mechanics/nodes/NodeGrid1D.h"
#include "nuto/mechanics/nodes/NodeDisplacements1D.h"

namespace NuTo
{
//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... class for grid nodes with displacements
class NodeGridDisplacements1D : public NodeGrid1D, public NodeDisplacements1D
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:

    //! @brief constructor
    NodeGridDisplacements1D(int rNodeId) : NodeGrid1D (rNodeId), NodeDisplacements1D()
    {}

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeGrid1D)
           & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeDisplacements1D);
    }
#endif  // ENABLE_SERIALIZATION

    //! @brief sets the global dofs
    //! @param rDOF current maximum DOF, this variable is increased within the routine
    virtual void SetGlobalDofs(int& rDOF)
    {
        NodeGrid1D::SetGlobalDofs(rDOF);
        NodeDisplacements1D::SetGlobalDofs(rDOF);
    }

    //! @brief write dof values to the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void SetGlobalDofValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
    {
        NodeGrid1D::SetGlobalDofValues(rActiveDofValues, rDependentDofValues);
        NodeDisplacements1D::SetGlobalDofValues(rActiveDofValues, rDependentDofValues);
    }

    //! @brief extract dof values from the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void GetGlobalDofValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
    {
        NodeGrid1D::GetGlobalDofValues(rActiveDofValues, rDependentDofValues);
        NodeDisplacements1D::GetGlobalDofValues(rActiveDofValues, rDependentDofValues);
    }

    //! @brief renumber the global dofs according to predefined ordering
    //! @param rMappingInitialToNewOrdering ... mapping from initial ordering to the new ordering
    virtual void RenumberGlobalDofs(std::vector<int>& rMappingInitialToNewOrdering)
    {
        NodeGrid1D::RenumberGlobalDofs(rMappingInitialToNewOrdering);
        NodeDisplacements1D::RenumberGlobalDofs(rMappingInitialToNewOrdering);
    }

    //! @brief returns the type of the node
    //! @return type
    virtual std::string GetNodeTypeStr()const
    {
    	return std::string("NodeGridDisplacements1D");
    }
};
}

#endif //NodeGridDisplacements_1d_H
