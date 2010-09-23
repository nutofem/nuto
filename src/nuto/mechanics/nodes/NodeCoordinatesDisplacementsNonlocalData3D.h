// $Id: $
#ifndef NodeCoordinatesDisplacementsNonlocalData3D_H
#define NodeCoordinatesDisplacementsNonlocalData3D_H

#include "nuto/mechanics/nodes/NodeCoordinatesDisplacements3D.h"
#include "nuto/mechanics/nodes/NodeNonlocalDataBase.h"

namespace NuTo
{
//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... class for nodes having coordinates and displacements
class NodeCoordinatesDisplacementsNonlocalData3D : public NodeCoordinatesDisplacements3D, public NodeNonlocalDataBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:

    //! @brief constructor
    NodeCoordinatesDisplacementsNonlocalData3D() : NodeCoordinatesDisplacements3D(), NodeNonlocalDataBase()
    {}

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeCoordinatesDisplacements3D)
           & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeNonlocalDataBase);
    }
#endif  // ENABLE_SERIALIZATION

};
}

#endif //NodeCoordinatesDisplacementsNonlocalData2D_H
