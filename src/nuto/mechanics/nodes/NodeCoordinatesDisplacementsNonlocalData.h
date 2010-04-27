// $Id: $
#ifndef NodeCoordinatesDisplacementsNonlocalData_H
#define NodeCoordinatesDisplacementsNonlocalData_H

#include "nuto/mechanics/nodes/NodeCoordinatesDisplacements.h"
#include "nuto/mechanics/nodes/NodeNonlocalDataBase.h"

namespace NuTo
{
//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... class for nodes having coordinates and displacements
template<int NUMCOORDINATES, int NUMDISPLACEMENTS>
class NodeCoordinatesDisplacementsNonlocalData : public  NodeCoordinatesDisplacements<NUMCOORDINATES,NUMDISPLACEMENTS>, public NodeNonlocalDataBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:

    //! @brief constructor
    NodeCoordinatesDisplacementsNonlocalData() : NodeCoordinatesDisplacements<NUMCOORDINATES,NUMDISPLACEMENTS> (), NodeNonlocalDataBase()
    {}

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeCoordinatesDisplacements)
           & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeNonlocalDataBase);
    }
#endif  // ENABLE_SERIALIZATION

};
}

#endif //NodeCoordinatesDisplacementsNonlocalData_H
