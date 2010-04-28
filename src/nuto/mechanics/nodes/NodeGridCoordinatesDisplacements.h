// $Id: $
#ifndef NodeGridCoordinatesDisplacements_H
#define NodeGridCoordinatesDisplacements_H

#include "nuto/mechanics/nodes/NodeGridCoordinates.h"
#include "nuto/mechanics/nodes/NodeDisplacements.h"

namespace NuTo
{
//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... class for nodes having coordinates and displacements
template<int NUMDISPLACEMENTS>
class NodeGridCoordinatesDisplacements : public  NodeGridCoordinates, public NodeDisplacements<NUMDISPLACEMENTS>
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:

    //! @brief constructor
    NodeGridCoordinatesDisplacements() : NodeGridCoordinates (), NodeDisplacements<NUMDISPLACEMENTS>()
    {}

   //! @brief constructor
    NodeGridCoordinatesDisplacements(unsigned int rNodeID) : NodeGridCoordinates (rNodeID), NodeDisplacements<NUMDISPLACEMENTS>()
    {}

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeGridCoordinates)
        & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeDisplacements);
    }
#endif  // ENABLE_SERIALIZATION

    //! @brief sets the global dofs
    //! @param rDOF current maximum DOF, this variable is increased within the routine
    virtual void SetGlobalDofs(int& rDOF)
    {
        NodeGridCoordinates::SetGlobalDofs(rDOF);
        NodeDisplacements<NUMDISPLACEMENTS>::SetGlobalDofs(rDOF);
    }

    //! @brief write dof values to the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void SetGlobalDofValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
    {
        NodeGridCoordinates::SetGlobalDofValues(rActiveDofValues, rDependentDofValues);
        NodeDisplacements<NUMDISPLACEMENTS>::SetGlobalDofValues(rActiveDofValues, rDependentDofValues);
    }

    //! @brief extract dof values from the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void GetGlobalDofValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
    {
        NodeGridCoordinates::GetGlobalDofValues(rActiveDofValues, rDependentDofValues);
        NodeDisplacements<NUMDISPLACEMENTS>::GetGlobalDofValues(rActiveDofValues, rDependentDofValues);
    }

    //! @brief renumber the global dofs according to predefined ordering
    //! @param rMappingInitialToNewOrdering ... mapping from initial ordering to the new ordering
    virtual void RenumberGlobalDofs(std::vector<int>& rMappingInitialToNewOrdering)
    {
        NodeGridCoordinates::RenumberGlobalDofs(rMappingInitialToNewOrdering);
        NodeDisplacements<NUMDISPLACEMENTS>::RenumberGlobalDofs(rMappingInitialToNewOrdering);
    }
};
}

#endif //NodeGridCoordinatesDisplacements_H
