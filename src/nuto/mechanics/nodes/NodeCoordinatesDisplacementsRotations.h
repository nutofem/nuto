// $Id: $
#ifndef NodeCoordinatesDisplacementsRotations_H
#define NodeCoordinatesDisplacementsRotations_H

#include "nuto/mechanics/nodes/NodeCoordinates.h"
#include "nuto/mechanics/nodes/NodeDisplacements.h"
#include "nuto/mechanics/nodes/NodeRotations.h"

namespace NuTo
{
//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... class for nodes having coordinates and displacements
template<int NUMCOORDINATES, int NUMDISPLACEMENTS, int NUMROTATIONS>
class NodeCoordinatesDisplacementsRotations : public  NodeCoordinates<NUMCOORDINATES>,
        public NodeDisplacements<NUMDISPLACEMENTS>,
        public NodeRotations<NUMROTATIONS>
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:

    //! @brief constructor
    NodeCoordinatesDisplacementsRotations() : NodeCoordinates<NUMCOORDINATES> (),
            NodeDisplacements<NUMDISPLACEMENTS>(),
            NodeRotations<NUMROTATIONS>()
    {}

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeCoordinates)
        & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeDisplacements)
        & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeRotations);
    }
#endif  // ENABLE_SERIALIZATION
    //! @brief sets the global dofs
    //! @param rDOF current maximum DOF, this variable is increased within the routine
    virtual void SetGlobalDofs(int& rDOF)
    {
        NodeCoordinates<NUMCOORDINATES>::SetGlobalDofs(rDOF);
        NodeDisplacements<NUMDISPLACEMENTS>::SetGlobalDofs(rDOF);
        NodeRotations<NUMROTATIONS>::SetGlobalDofs(rDOF);
    }

    //! @brief write dof values to the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void SetGlobalDofValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
    {
        NodeCoordinates<NUMCOORDINATES>::SetGlobalDofValues(rActiveDofValues, rDependentDofValues);
        NodeDisplacements<NUMDISPLACEMENTS>::SetGlobalDofValues(rActiveDofValues, rDependentDofValues);
        NodeRotations<NUMROTATIONS>::SetGlobalDofValues(rActiveDofValues, rDependentDofValues);
    }

    //! @brief extract dof values from the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void GetGlobalDofValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
    {
        NodeCoordinates<NUMCOORDINATES>::GetGlobalDofValues(rActiveDofValues, rDependentDofValues);
        NodeDisplacements<NUMDISPLACEMENTS>::GetGlobalDofValues(rActiveDofValues, rDependentDofValues);
        NodeRotations<NUMROTATIONS>::GetGlobalDofValues(rActiveDofValues, rDependentDofValues);
    }

    //! @brief renumber the global dofs according to predefined ordering
    //! @param rMappingInitialToNewOrdering ... mapping from initial ordering to the new ordering
    virtual void RenumberGlobalDofs(std::vector<int>& rMappingInitialToNewOrdering)
    {
        NodeCoordinates<NUMCOORDINATES>::RenumberGlobalDofs(rMappingInitialToNewOrdering);
        NodeDisplacements<NUMDISPLACEMENTS>::RenumberGlobalDofs(rMappingInitialToNewOrdering);
        NodeRotations<NUMROTATIONS>::RenumberGlobalDofs(rMappingInitialToNewOrdering);
    }
};
}

#endif //NodeCoordinatesDisplacementsRotations_H
