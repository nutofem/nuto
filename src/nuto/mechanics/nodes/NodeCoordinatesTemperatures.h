// $Id: $
#ifndef NodeCoordinatesTemperatures_H
#define NodeCoordinatesTemperatures_H

#include "nuto/mechanics/nodes/NodeCoordinates.h"
#include "nuto/mechanics/nodes/NodeTemperatures.h"

namespace NuTo
{
//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... class for nodes having coordinates and Temperatures
template<int NUMCOORDINATES, int NUMTEMPERATURES>
class NodeCoordinatesTemperatures : public  NodeCoordinates<NUMCOORDINATES>, public NodeTemperatures<NUMTEMPERATURES>
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:

    //! @brief constructor
    NodeCoordinatesTemperatures() : NodeCoordinates<NUMCOORDINATES> (), NodeTemperatures<NUMTEMPERATURES>()
    {}

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeCoordinates)
        & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeTemperatures);
    }
#endif  // ENABLE_SERIALIZATION
    //! @brief sets the global dofs
    //! @param rDOF current maximum DOF, this variable is increased within the routine
    virtual void SetGlobalDofs(int& rDOF)
    {
        NodeCoordinates<NUMCOORDINATES>::SetGlobalDofs(rDOF);
        NodeTemperatures<NUMTEMPERATURES>::SetGlobalDofs(rDOF);
    }

    //! @brief write dof values to the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void SetGlobalDofValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
    {
        NodeCoordinates<NUMCOORDINATES>::SetGlobalDofValues(rActiveDofValues, rDependentDofValues);
        NodeTemperatures<NUMTEMPERATURES>::SetGlobalDofValues(rActiveDofValues, rDependentDofValues);
    }

    //! @brief extract dof values from the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void GetGlobalDofValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
    {
        NodeCoordinates<NUMCOORDINATES>::GetGlobalDofValues(rActiveDofValues, rDependentDofValues);
        NodeTemperatures<NUMTEMPERATURES>::GetGlobalDofValues(rActiveDofValues, rDependentDofValues);
    }

    //! @brief renumber the global dofs according to predefined ordering
    //! @param rMappingInitialToNewOrdering ... mapping from initial ordering to the new ordering
    virtual void RenumberGlobalDofs(std::vector<int>& rMappingInitialToNewOrdering)
    {
        NodeCoordinates<NUMCOORDINATES>::RenumberGlobalDofs(rMappingInitialToNewOrdering);
        NodeTemperatures<NUMTEMPERATURES>::RenumberGlobalDofs(rMappingInitialToNewOrdering);
    }
};
}

#endif //NodeCoordinatesTemperatures_H
